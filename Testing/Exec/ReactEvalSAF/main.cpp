#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <PltFileManager.H>

#ifdef AMREX_USE_GPU
#include <AMReX_SUNMemory.H>
#endif

#include "mechanism.H"
#include <GPU_misc.H>

#include <PelePhysics.H>
#include <ReactorBase.H>

namespace {
const std::string level_prefix{"Level_"};
}

void
GotoNextLine(std::istream& is)
{
  constexpr std::streamsize bl_ignore_max{100000};
  is.ignore(bl_ignore_max, '\n');
}

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  {

    amrex::Real strt_time = amrex::ParallelDescriptor::second();

    BL_PROFILE_VAR("main::main()", pmain);

    amrex::ParmParse pp;
    std::string fuel_name;
    pp.get("fuel_name", fuel_name);

    std::string chem_integrator = "";
    pp.get("chem_integrator", chem_integrator);

    std::string pltfile;
    bool do_plt = false;
    if (pp.countval("plotfile") > 0) {
      pp.get("plotfile", pltfile);
      do_plt = true;
    }

    /* initial data */
    // Can either generate single level data
    // or read in data from a plotfile like file
    int initFromPlt = 0;
    int enforce_forcing_zero=0;
    pp.query("initFromFile", initFromPlt);
    std::string pltfileIn = "";
    if (initFromPlt) {
      pp.query("initFile", pltfileIn);
    }

    bool subcycledt = true;
    pp.query("subcycle_dt", subcycledt);

    /* react() function version */
    // 1 -> Array4 version of react()  (Default)
    // 2 -> 1d raw pointer version of react()
    std::string reactFormat = "Array4";
    int reactFunc;
    pp.query("reactFormat", reactFormat);
    if (reactFormat == "Array4") {
      reactFunc = 1;
    } else if (reactFormat == "1dArray") {
      reactFunc = 2;
    } else {
      amrex::Abort(" --> reactFormat can only be 'Array4' or '1dArray' !");
    }

    /* ODE inputs */
    amrex::ParmParse ppode("ode");
    int ode_ncells = 1;
    ppode.query("ode_ncells", ode_ncells); // number of cells to integrate per
                                           // call, used only if reactFunc = 2

    amrex::Real dt = 1.e-5;
    ppode.query("dt", dt);

    ppode.query("enforce_forcing_zero",enforce_forcing_zero);

    int ndt = 1;
    ppode.query("ndt", ndt); // number of solver calls per dt

    int ode_iE = -1;
    ppode.query(
      "reactor_type",
      ode_iE); // RHS type, 1: e (PeleC), !1: h (PeleLM)  <------ FIXME!

    amrex::Real rtol = 1e-10;
    ppode.query("rtol", rtol);

    amrex::Real atol = 1e-10;
    ppode.query("atol", atol);

    int use_typ_vals = 0;
    ppode.query("use_typ_vals", use_typ_vals);

    amrex::Print() << "ODE solver: " << chem_integrator << std::endl;
    amrex::Print() << "Type of reactor: "
                   << (ode_iE == 1 ? "e (PeleC)" : "h (PeleLM)")
                   << std::endl; // <---- FIXME
    amrex::Print() << "Fuel: " << fuel_name << ", Oxy: O2" << std::endl;

    /* Mixture info */
    int fuel_idx = -1;
    if (fuel_name == "H2") {
      fuel_idx = H2_ID;
#ifdef CH4_ID
    } else if (fuel_name == "CH4") {
      fuel_idx = CH4_ID;
#endif
#ifdef NC12H26_ID
    } else if (fuel_name == "NC12H26") {
      fuel_idx = NC12H26_ID;
#endif
#ifdef IC8H18_ID
    } else if (fuel_name == "IC8H18") {
      fuel_idx = IC8H18_ID;
#endif
    }

    /* Initialize reactor object inside OMP region, including tolerances */
    BL_PROFILE_VAR("main::reactor_info()", reactInfo);
    std::unique_ptr<pele::physics::reactions::ReactorBase> reactor =
      pele::physics::reactions::ReactorBase::create(chem_integrator);
    reactor->init(ode_iE, ode_ncells);
    BL_PROFILE_VAR_STOP(reactInfo);

    // -----------------------------------------------------------------------------
    // Initialize geom/data
    // When initFromPlt = 0, default is single level data provided by
    // initialize_data with a -1:1 unit length realbox in each dir and cell
    // count provided by the user otherwise number of levels, geom and all are
    // read in from the pltfileIn
    // -----------------------------------------------------------------------------
    int finest_level = 0;
    amrex::Vector<amrex::Geometry> geoms;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmaps;

    BL_PROFILE_VAR("main::geometry_setup", GeomSetup);
    if (initFromPlt == 0) {
      // -----------------------------------------------------------------------------
      // Resize vectors
      // -----------------------------------------------------------------------------
      geoms.resize(finest_level + 1);
      grids.resize(finest_level + 1);
      dmaps.resize(finest_level + 1);

      // -----------------------------------------------------------------------------
      // Define the geometry
      // -----------------------------------------------------------------------------
      std::array<int, 3> ncells = {AMREX_D_DECL(1, 1, 1)};
      if (pp.countval("ncells") == 1) {
        pp.get("ncells", ncells[0]);
        ncells = {AMREX_D_DECL(ncells[0], 1, 1)};
      } else if (pp.countval("ncells") >= AMREX_SPACEDIM) {
        amrex::Vector<int> nc(AMREX_SPACEDIM);
        pp.getarr("ncells", nc, 0, AMREX_SPACEDIM);
        ncells = {AMREX_D_DECL(nc[0], nc[1], nc[2])};
      } else {
        amrex::Abort("ncells has to have length 1 or spacedim");
      }

      amrex::Box domain(
        amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
        amrex::IntVect(
          AMREX_D_DECL(ncells[0] - 1, ncells[1] - 1, ncells[2] - 1)));

      amrex::RealBox real_box(
        {AMREX_D_DECL(-1.0, -1.0, -1.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});

      int coord = 0;

      amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};

      geoms[0] = amrex::Geometry(domain, real_box, coord, is_periodic);

      amrex::Print() << "Integrating " << domain.numPts()
                     << " cells for: " << dt << " seconds with " << ndt
                     << " substeps \n";

      // -----------------------------------------------------------------------------
      // Define amrex::BoxArray / Dmap
      // -----------------------------------------------------------------------------
      int max_grid_size = 16;
      pp.query("max_grid_size", max_grid_size);
      grids[0] = amrex::BoxArray(domain);
      grids[0].maxSize(max_grid_size);
      dmaps[0] = amrex::DistributionMapping(
        grids[0], amrex::ParallelDescriptor::NProcs());
    } else {
      // Use PltFileManager, only metadata loaded here
      pele::physics::pltfilemanager::PltFileManager pltData(pltfileIn);

      if (ode_iE == 1) {
        amrex::Abort(
          "The option to read in chemical data is currently only available "
          "with PeleLM data and requires ode_iE=2");
      }

      //--------------------------------------------
      // General info
      //--------------------------------------------
      finest_level = pltData.getNlev() - 1;
      //dt = pltData.getTime(); // In LMeX we stored the dt in the pltfileIn time
      amrex::Print()<<"\ndt read from plt file = "<<dt;

      // -----------------------------------------------------------------------------
      // Resize vectors
      // -----------------------------------------------------------------------------
      geoms.resize(finest_level + 1);
      grids.resize(finest_level + 1);
      dmaps.resize(finest_level + 1);

      // -----------------------------------------------------------------------------
      // Define geoms, read amrex::BoxArray and define dmap
      // -----------------------------------------------------------------------------
      geoms[0] = pltData.getGeom(0);
      for (int lev = 1; lev <= finest_level; ++lev) {
        geoms[lev] = amrex::refine(geoms[lev - 1], 2); // Assumes ref_ratio = 2
      }

      for (int lev = 0; lev <= finest_level; ++lev) {
        // Set vector entries
        grids[lev] = pltData.getGrid(lev);
        dmaps[lev] = amrex::DistributionMapping(
          grids[lev], amrex::ParallelDescriptor::NProcs());
      }

      for (int lev = 0; lev <= finest_level; ++lev) {
        if (subcycledt) {
          amrex::Print() << "  Level " << lev << " integrating "
                         << grids[lev].numPts() << " cells on "
                         << grids[lev].size() << " boxes for "
                         << dt / std::pow(2, lev) << " seconds \n";
        } else {
          amrex::Print() << "  Level " << lev << " integrating "
                         << grids[lev].numPts() << " cells on "
                         << grids[lev].size() << " boxes for " << dt
                         << " seconds \n";
        }
      }
    }
    BL_PROFILE_VAR_STOP(GeomSetup);

    // -----------------------------------------------------------------------------
    // Create MFs and generate initial data
    // -----------------------------------------------------------------------------
    BL_PROFILE_VAR("main::initialize_data()", InitData);
    int num_grow = 0;
    amrex::Vector<amrex::MultiFab> mf(finest_level + 1);
    amrex::Vector<amrex::MultiFab> rY_source_ext(finest_level + 1);
    amrex::Vector<amrex::MultiFab> mfE(finest_level + 1);
    amrex::Vector<amrex::MultiFab> rY_source_energy_ext(finest_level + 1);
    amrex::Vector<amrex::MultiFab> fctCount(finest_level + 1);
    amrex::Vector<amrex::iMultiFab> dummyMask(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
      mf[lev].define(grids[lev], dmaps[lev], NUM_SPECIES + 1, num_grow);
      rY_source_ext[lev].define(grids[lev], dmaps[lev], NUM_SPECIES, num_grow);
      mfE[lev].define(grids[lev], dmaps[lev], 1, num_grow);
      rY_source_energy_ext[lev].define(grids[lev], dmaps[lev], 1, num_grow);
      fctCount[lev].define(grids[lev], dmaps[lev], 1, num_grow);
      dummyMask[lev].define(grids[lev], dmaps[lev], 1, num_grow);
      dummyMask[lev].setVal(1);
    }

    amrex::FabArrayBase::mfiter_tile_size =
      amrex::IntVect(AMREX_D_DECL(1024, 1024, 1024));

    // -----------------------------------------------------------------------------
    // Initialize data
    // -----------------------------------------------------------------------------
    if (initFromPlt == 0) {
      for (int lev = 0; lev <= finest_level; ++lev) {
        const auto geomdata = geoms[lev].data();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(mf[lev], amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {

          const amrex::Box& box = mfi.tilebox();

          amrex::Array4<amrex::Real> const& rY_a = mf[lev].array(mfi);
          amrex::Array4<amrex::Real> const& rYs_a =
            rY_source_ext[lev].array(mfi);
          amrex::Array4<amrex::Real> const& E_a = mfE[lev].array(mfi);
          amrex::Array4<amrex::Real> const& rE_a =
            rY_source_energy_ext[lev].array(mfi);

          amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              initialize_data(
                i, j, k, fuel_idx, ode_iE, rY_a, rYs_a, E_a, rE_a, geomdata);
            });
        }
      }
    } else {
      pele::physics::pltfilemanager::PltFileManager pltData(pltfileIn);
      amrex::Vector<std::string> plt_vars = pltData.getVariableList();

      // Figure out the indices of reauired entries
      // Assuming PeleLMeX evaluate mode was used: rhoY + rhoH + T + F_rhoY +
      // F_rhoH
      amrex::Vector<std::string> spec_names;
      pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
        spec_names);
      const std::string firstSpecVar = "rho.Y(" + spec_names[0] + ")";
      const std::string firstSpecForceVar = "F_rhoY(" + spec_names[0] + ")";
      int idFirstSpec = -1, idFirstSpecForce = -1, idrhoH = -1, idT = -1,
          idrhoHForce = -1;
      for (int i = 0; i < plt_vars.size(); ++i) {
    	  amrex::Print()<<"\nPlt_vars = "<<plt_vars[i];
        if (plt_vars[i] == firstSpecVar) {
          idFirstSpec = i;
        } else if (plt_vars[i] == firstSpecForceVar) {
          idFirstSpecForce = i;
        } else if (plt_vars[i] == "rhoh") {
          idrhoH = i;
        } else if (plt_vars[i] == "temp") {
          idT = i;
        } else if (plt_vars[i] == "F_rhoH") {
          idrhoHForce = i;
        }
      }
      amrex::Print()<<"\n Values = "<<idFirstSpec<<" "<<idFirstSpecForce<<" "<<idrhoH<<" "<<idT<<" "<<idrhoHForce;
      if (idFirstSpec * idFirstSpecForce * idrhoH * idT * idrhoHForce < 0) {

        amrex::Abort("All of 'rhoY(<spec>)', 'F_rhoY(<spec>)', 'rhoH', 'Temp' "
                     "and 'F_rhoH' must be in the pltfile !");
      }

      // Load the data from pltfile
      for (int lev = 0; lev <= finest_level; ++lev) {
        // rhoYs
        pltData.fillPatchFromPlt(
          lev, geoms[lev], idFirstSpec, 0, NUM_SPECIES, mf[lev]);

        // Temperature
        pltData.fillPatchFromPlt(lev, geoms[lev], idT, NUM_SPECIES, 1, mf[lev]);

        // RhoH
        pltData.fillPatchFromPlt(lev, geoms[lev], idrhoH, 0, 1, mfE[lev]);

        // Force rhoYs
        pltData.fillPatchFromPlt(
          lev, geoms[lev], idFirstSpecForce, 0, NUM_SPECIES,
          rY_source_ext[lev]);

        // Force rhoH
        pltData.fillPatchFromPlt(
          lev, geoms[lev], idrhoHForce, 0, 1, rY_source_energy_ext[lev]);

        // Convert from MKS -> CGS since we have PeleLMeX data
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(mf[lev], amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
          const amrex::Box& box = mfi.tilebox();

          auto const& rhoYs = mf[lev].array(mfi);
          auto const& temp = mf[lev].array(mfi, NUM_SPECIES);
          auto const& rhoH = mfE[lev].array(mfi);
          auto const& force_Y = rY_source_ext[lev].array(mfi);
          auto const& force_E = rY_source_energy_ext[lev].array(mfi);
          amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              for (int n = 0; n < NUM_SPECIES; n++) {
                //Sreejith;s change
            	 amrex::Real problematic_y[] =  {0,	0.00000001906452056,	0.0000004527277204,	0.000003901857124,	0.000000005778085447,	0.0000001890299964,	0.00009508136884,	0.000000000550059528,	0.00002160084221,	0,	0,	0,	0.00001562552901,	0.0002242665362,	0,	0,	0,	0,	1.038992968E-11,	3.741023963E-12,	1.47186468E-12,	1.197974087E-13,	0,	0,	0,	0,	3.526834578E-15,	0,	0,	0,	0,	0,	0,	0,	0.000752485258,	-1310669.84};
            	//amrex::Real problematic_y[] = {-2.212600914E-14,	0.00000001923291707,	0.0000004545332138,	0.000003906022951,	0.000000005812144857,	0.0000001902028235,	0.00009504791236,	5.521241896E-10,	0.0000216376225,	8.431245768E-12,	4.40519813E-12,	-4.362354034E-12,	0.00001569572293,	0.0002241411391,	3.339259781E-11,	1.259068579E-11,	-2.667453863E-13,	-2.938625592E-13,	1.896291021E-11,	4.822275113E-12,	1.189467809E-12,	2.057710371E-12,	-2.342593708E-12,	9.716986658E-16,	2.463168709E-13,	2.755314714E-13,	3.414364215E-13,	7.482578953E-13,	4.672226628E-13,	-2.568778166E-14,	3.506010798E-13,	3.764985622E-14,	-7.362617213E-18,	-4.905050178E-18,	0.0007526698623,	2456.330338};
                amrex::Real problematic_ysource[] = {-0.0002597907082,	0.001622441092,	0.03082477837,	0.001612022016,	0.00007165398376,	0.011534525,	-0.2662386564,	0.00002241043973,	0.3293143519,	0.0001157438953,	0.00003746859982,	-0.0001083773136,	0.5853110608,	-1.03563453,	0.0001140656165,	0.0002988066165,	-0.000007944280146,	-0.000005147492599,	0.00008608740849,	0.0002040468676,	0.0001050412546,	0.0000299641377,	-0.0000915634789,	-3.955449377E-11,	0.000002302668235,	0.000002337397766,	0.000003129288029,	0.000006240488129,	0.000004209584895,	0.0000003151098105,	0.000002892875026,	-0.000000483061241,	-0.000002215056785,	-0.000001597039545,	1.415814843};
                rhoYs(i, j, k, n) = problematic_y[n]* 1.0;
                force_Y(i, j, k, n) = problematic_ysource[n]* 1.0;


            	/*rhoYs(i, j, k, n) *= 0.001;
            	force_Y(i, j, k, n) *= 0.001;*/
                if(enforce_forcing_zero==1)
                	{
                	force_Y(i, j, k, n) = 0.0;		//Sreejith. Enforcing forcing to zero
                	}
              }
              amrex::Real problematic_temp = 2456.330338;
              temp(i,j,k) = problematic_temp;

              amrex::Real problematic_energy = -1310669.84 ;
              amrex::Real problematic_energy_source = 7.16431189500000E+10;
              rhoH(i, j, k) = problematic_energy*1.0;
              force_E(i, j, k) = problematic_energy_source* 1.0;

              /*rhoH(i, j, k) *= 10.0;
              force_E(i, j, k) *= 10.0;*/

              if(enforce_forcing_zero==1)
              {
            	  force_E(i, j, k) = 0.0;		//Sreejith. Enforcing forcing to zero
               }
            });
        }
      }
    }
    BL_PROFILE_VAR_STOP(InitData);

    // Print initial data if needed
    BL_PROFILE_VAR_NS("PlotFile", PlotFile);
    if (do_plt) {
      BL_PROFILE_VAR_START(PlotFile);
      std::string outfile = amrex::Concatenate(pltfile, 0);
      amrex::Vector<int> isteps(finest_level + 1, 0);
      amrex::Vector<amrex::IntVect> refRatios(
        finest_level, {AMREX_D_DECL(2, 2, 2)});
      amrex::Vector<std::string> plt_VarsName;
      for (int k = 0; k < NUM_SPECIES; ++k) {
        plt_VarsName.push_back("SPEC" + std::to_string(k));
      }
      plt_VarsName.push_back("TEMP");

      amrex::WriteMultiLevelPlotfile(
        outfile, finest_level + 1, GetVecOfConstPtrs(mf), plt_VarsName, geoms,
        0.0, isteps, refRatios);
      BL_PROFILE_VAR_STOP(PlotFile);
    }

    // -----------------------------------------------------------------------------
    // Set typical values
    // -----------------------------------------------------------------------------
    if (
      chem_integrator == "ReactorCvode" || chem_integrator == "ReactorArkode") {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
      {
        if (use_typ_vals) {
          amrex::Vector<amrex::Real> typ_vals(NUM_SPECIES + 1, 1.0e-10);
          if (ppode.contains("typ_vals")) {
            amrex::Print()
              << "Using user-defined typical values for the absolute "
                 "tolerances of the ode solver.\n";
            ppode.getarr("typ_vals", typ_vals, 0, NUM_SPECIES + 1);
            for (int i = 0; i < NUM_SPECIES; ++i) {
              typ_vals[i] = std::max(typ_vals[i], 1.e-10);
            }
          } else {
            amrex::Print()
              << "Using typical values from the initial data for the "
                 "absolute tolerances of the ode solver.\n";
            for (int lev = 0; lev <= finest_level; ++lev) {
              /*
                Pretty sure TypVal should be rhoYs in CGS, but keep that around
                just in case. amrex::MultiFab
                massFrac(grids[lev],dmaps[lev],NUM_SPECIES,0); #ifdef
                AMREX_USE_OMP #pragma omp parallel if
                (amrex::Gpu::notInLaunchRegion()) #endif for (amrex::MFIter
                mfi(mf[lev],amrex::TilingIfNotGPU()); mfi.isValid();
                ++mfi) { const amrex::Box& box = mfi.tilebox(); auto const&
                rhoYs = mf[lev].const_array(mfi); auto const& Ys    =
                massFrac.array(mfi); amrex::ParallelFor(box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                amrex::Real rho = 0.0;
                for (int n = 0; n < NUM_SPECIES; n++) {
                rho += rhoYs(i,j,k,n);
                }
                amrex::Real rhoinv = 1.0/rho;
                for (int n = 0; n < NUM_SPECIES; n++) {
                Ys(i,j,k,n) = rhoYs(i,j,k,n) * rhoinv;
                }
                });
                }
              */
              for (int i = 0; i < NUM_SPECIES; ++i) {
                typ_vals[i] = std::max(typ_vals[i], mf[lev].max(i));
                amrex::Print()<<"\nWhile calculating typical values for i = "<<i<<" max = "<<mf[lev].max(i)<<" min = "<<mf[lev].min(i);
              }
              typ_vals[NUM_SPECIES] =
                std::max(typ_vals[NUM_SPECIES], mf[lev].max(NUM_SPECIES));
            }
          }
          reactor->set_typ_vals_ode(typ_vals);
        }
      }
    }
    amrex::Print() << " \n STARTING THE ADVANCE \n";

    /* REACT */
    BL_PROFILE_VAR_NS("React", ReactInLoop);
    BL_PROFILE_VAR_NS("Allocs", Allocs);
    BL_PROFILE_VAR_NS("Flatten", mainflatten);

    for (int lev = 0; lev <= finest_level; ++lev) {
      amrex::Real lvl_strt = amrex::ParallelDescriptor::second();
      BL_PROFILE_VAR("Advance_Level" + std::to_string(lev), Advance);
#ifdef AMREX_USE_OMP
      const auto tiling = amrex::MFItInfo().SetDynamic(true);
#pragma omp parallel
#else
      const bool tiling = amrex::TilingIfNotGPU();
#endif
      for (amrex::MFIter mfi(mf[lev], tiling); mfi.isValid(); ++mfi) {

        int omp_thread = 0;
#ifdef AMREX_USE_OMP
        omp_thread = omp_get_thread_num();
#endif
        const amrex::Box& box = mfi.tilebox();
        const int nc = box.numPts();
        int extra_cells = 0;

        auto const& rhoY = mf[lev].array(mfi);
        auto const& T = mf[lev].array(mfi, NUM_SPECIES);
        auto const& rhoE = mfE[lev].array(mfi);
        auto const& frcExt = rY_source_ext[lev].array(mfi);
        auto const& frcEExt = rY_source_energy_ext[lev].array(mfi);
        auto const& fc = fctCount[lev].array(mfi);
        auto const& mask = dummyMask[lev].array(mfi);
        // -------------------------------------------------------------
        // Integration with Array4 react function
        if (reactFunc == 1) {
          amrex::Real time = 0.0;
          amrex::Real dt_lev = (subcycledt) ? dt / std::pow(2, lev) : dt;
          amrex::Real dt_incr = dt_lev / ndt;
          int tmp_fc;
          //if (omp_thread == 0) {
            amrex::Print() << "  [" << lev << "]"
                           << " integrating " << nc << " cells using dt = "<<dt_incr<<" \n";
          //}
          /* Solve */
          BL_PROFILE_VAR_START(ReactInLoop);
          for (int ii = 0; ii < ndt; ++ii) {
            tmp_fc = reactor->react(
              box, rhoY, frcExt, T, rhoE, frcEExt, fc, mask, dt_incr, time
#ifdef AMREX_USE_GPU
              ,
              amrex::Gpu::gpuStream()
#endif
            );
            dt_incr = dt_lev / ndt;
            amrex::Gpu::Device::streamSynchronize();
          }
          BL_PROFILE_VAR_STOP(ReactInLoop);

          // -------------------------------------------------------------
          // Integration with 1dArray raw pointer react function
        } else if (reactFunc == 2) {

          // On GPU, integrate the entirely box at once
          // othewise use the user-input ode_ncells
#ifdef AMREX_USE_GPU
          ode_ncells = nc;
#endif
          extra_cells = nc - (nc / ode_ncells) * ode_ncells;
          amrex::Print()<<"\n Extra cells = "<<extra_cells;
          if (omp_thread == 0) {
            amrex::Print() << " Integrating " << nc << " cells with a "
                           << ode_ncells << " ode cell buffer ";
            amrex::Print() << "(" << extra_cells << " extra cells) \n";
          }

          BL_PROFILE_VAR_START(Allocs);
          int nCells = nc + extra_cells;

#ifdef AMREX_USE_GPU
          auto tmp_vect_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
            nCells * (NUM_SPECIES + 1) * sizeof(amrex::Real));
          auto tmp_src_vect_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
            nCells * NUM_SPECIES * sizeof(amrex::Real));
          auto tmp_vect_energy_d =
            (amrex::Real*)amrex::The_Device_Arena()->alloc(
              nCells * sizeof(amrex::Real));
          auto tmp_src_vect_energy_d =
            (amrex::Real*)amrex::The_Device_Arena()->alloc(
              nCells * sizeof(amrex::Real));
          auto tmp_fc_d = (long int*)amrex::The_Device_Arena()->alloc(
            nCells * sizeof(amrex::Real));
          auto tmp_mask_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
            nCells * sizeof(amrex::Real));
#endif

          auto tmp_vect = new amrex::Real[nCells * (NUM_SPECIES + 1)];
          auto tmp_src_vect = new amrex::Real[nCells * NUM_SPECIES];
          auto tmp_vect_energy = new amrex::Real[nCells];
          auto tmp_src_vect_energy = new amrex::Real[nCells];
          auto tmp_fc = new long int[nCells];
          auto tmp_mask = new int[nCells];

          BL_PROFILE_VAR_STOP(Allocs);

          BL_PROFILE_VAR_START(mainflatten);
#ifndef AMREX_USE_GPU
          reactor->flatten(
            box, nCells, rhoY, frcExt, T, rhoE, frcEExt, tmp_vect, tmp_src_vect,
            tmp_vect_energy, tmp_src_vect_energy);

          //Forget the mf values. Use these values I got from running Pele on frontier

          for (int icell = nc; icell < nc + extra_cells; icell++) {
            for (int sp = 0; sp < NUM_SPECIES; sp++) {
              tmp_vect[icell * (NUM_SPECIES + 1) + sp] = rhoY(0, 0, 0, sp);
              tmp_src_vect[icell * NUM_SPECIES + sp] = frcExt(0, 0, 0, sp);
            }
            tmp_vect[icell * (NUM_SPECIES + 1) + NUM_SPECIES] = T(0, 0, 0);
            tmp_vect_energy[icell] = rhoE(0, 0, 0);
            tmp_src_vect_energy[icell] = frcEExt(0, 0, 0);
            tmp_mask[icell] = mask(0, 0, 0);
          }
#else
          reactor->flatten(
            box, nCells, rhoY, frcExt, T, rhoE, frcEExt, tmp_vect_d,
            tmp_src_vect_d, tmp_vect_energy_d, tmp_src_vect_energy_d);

          amrex::Gpu::copy(
            amrex::Gpu::deviceToHost, tmp_vect_d,
            tmp_vect_d + nCells * (NUM_SPECIES + 1), tmp_vect);
          amrex::Gpu::copy(
            amrex::Gpu::deviceToHost, tmp_src_vect_d,
            tmp_src_vect_d + nCells * NUM_SPECIES, tmp_src_vect);
          amrex::Gpu::copy(
            amrex::Gpu::deviceToHost, tmp_vect_energy_d,
            tmp_vect_energy_d + nCells, tmp_vect_energy);
          amrex::Gpu::copy(
            amrex::Gpu::deviceToHost, tmp_src_vect_energy_d,
            tmp_src_vect_energy_d + nCells, tmp_src_vect_energy);
#endif
          BL_PROFILE_VAR_STOP(mainflatten);

          /* Solve */
          BL_PROFILE_VAR_START(ReactInLoop);
          for (int i = 0; i < nCells; i += ode_ncells) {
        	  amrex::Print()<<"\n Solving for cell "<<i;
            tmp_fc[i] = 0;
            amrex::Real time = 0.0;
            amrex::Real dt_lev = (subcycledt) ? dt / std::pow(2, lev) : dt;
            amrex::Real dt_incr = dt_lev / ndt;
            amrex::Print()<<"\nInternal loop "<<dt_incr<<" "<<time;
            for (int ii = 0; ii < ndt; ++ii) {



              tmp_fc[i] += reactor->react(
                &tmp_vect[i * (NUM_SPECIES + 1)],
                &tmp_src_vect[i * NUM_SPECIES], &tmp_vect_energy[i],
                &tmp_src_vect_energy[i], dt_incr, time, ode_ncells
#ifdef AMREX_USE_GPU
                ,
                amrex::Gpu::gpuStream()
#endif
              );
              amrex::Print()<<"\n Solved for cell "<<i;

              dt_incr = dt_lev / ndt;
              for (int ic = i + 1; ic < i + ode_ncells; ++ic) {
                tmp_fc[ic] = tmp_fc[i];
              }
              amrex::Gpu::Device::streamSynchronize();
            }
          }
          BL_PROFILE_VAR_STOP(ReactInLoop);

          BL_PROFILE_VAR_START(mainflatten);
#ifndef AMREX_USE_GPU
          reactor->unflatten(
            box, nCells, rhoY, T, rhoE, frcEExt, fc, tmp_vect, tmp_vect_energy,
            tmp_fc, dt);
#else

          amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, tmp_vect,
            tmp_vect + nCells * (NUM_SPECIES + 1), tmp_vect_d);
          amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, tmp_src_vect,
            tmp_src_vect + nCells * NUM_SPECIES, tmp_src_vect_d);
          amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, tmp_vect_energy, tmp_vect_energy + nCells,
            tmp_vect_energy_d);
          amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, tmp_src_vect_energy,
            tmp_src_vect_energy + nCells, tmp_src_vect_energy_d);
          amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, tmp_fc, tmp_fc + nCells, tmp_fc_d);

          reactor->unflatten(
            box, nCells, rhoY, T, rhoE, frcEExt, fc, tmp_vect_d,
            tmp_vect_energy_d, tmp_fc_d, dt);
#endif
          BL_PROFILE_VAR_STOP(mainflatten);

          delete[] tmp_vect;
          delete[] tmp_src_vect;
          delete[] tmp_vect_energy;
          delete[] tmp_src_vect_energy;
          delete[] tmp_fc;
          delete[] tmp_mask;
#ifdef AMREX_USE_GPU
          amrex::The_Device_Arena()->free(tmp_vect_d);
          amrex::The_Device_Arena()->free(tmp_src_vect_d);
          amrex::The_Device_Arena()->free(tmp_vect_energy_d);
          amrex::The_Device_Arena()->free(tmp_src_vect_energy_d);
          amrex::The_Device_Arena()->free(tmp_fc_d);
          amrex::The_Device_Arena()->free(tmp_mask_d);
#endif
        }
      }
      BL_PROFILE_VAR_STOP(Advance);
      amrex::Real lvl_run_time = amrex::ParallelDescriptor::second() - lvl_strt;
      amrex::ParallelDescriptor::ReduceRealMax(
        lvl_run_time, amrex::ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << "   >> Level " << lev << " advance: " << lvl_run_time
                     << "\n";
    }

    // TODO multilevel max.
    {
      amrex::Vector<double> typ_vals(NUM_SPECIES + 1);
      amrex::Print() << "ode.typ_vals= ";
      for (int i = 0; i < NUM_SPECIES + 1; ++i) {
        amrex::Print() << std::max(1.e-10, mf[0].max(i)) << " ";
      }
      amrex::Print() << std::endl;
    }

    if (do_plt) {
      BL_PROFILE_VAR_START(PlotFile);
      std::string outfile = amrex::Concatenate(pltfile, 1);
      // TODO: add fct count to this output
      amrex::Vector<int> isteps(finest_level + 1, 0);
      amrex::Vector<amrex::IntVect> refRatios(
        finest_level, {AMREX_D_DECL(2, 2, 2)});
      amrex::Vector<std::string> plt_VarsName;
      for (int k = 0; k < NUM_SPECIES; ++k) {
        plt_VarsName.push_back("SPEC" + std::to_string(k));
      }
      plt_VarsName.push_back("TEMP");

      amrex::WriteMultiLevelPlotfile(
        outfile, finest_level + 1, GetVecOfConstPtrs(mf), plt_VarsName, geoms,
        0.0, isteps, refRatios);
      BL_PROFILE_VAR_STOP(PlotFile);
    }

    BL_PROFILE_VAR_STOP(pmain);

    amrex::Real run_time = amrex::ParallelDescriptor::second() - strt_time;
    amrex::ParallelDescriptor::ReduceRealMax(
      run_time, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << " \n >> React::main() " << run_time << "\n\n";
  }

  amrex::Finalize();

  return 0;
}
