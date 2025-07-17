#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include "Utilities.H"

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);



  // Test Unit Conversions
  namespace m2c = pele::physics::utilities::mks2cgs;
  namespace c2m = pele::physics::utilities::cgs2mks;

    // Define geometry
    amrex::Array<int, AMREX_SPACEDIM> npts{AMREX_D_DECL(1, 1, 1)};
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      npts[i] = 16;
    }

    amrex::Box domain(
      amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
      amrex::IntVect(AMREX_D_DECL(npts[0] - 1, npts[1] - 1, npts[2] - 1)));
    int max_size = 32;
    amrex::BoxArray ba(domain);
    ba.maxSize(max_size);
    amrex::RealBox real_box(
      {AMREX_D_DECL(-1.0, -1.0, -1.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});

    int coord = 0;
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
    amrex::Geometry geom(domain, real_box, coord, is_periodic);
    amrex::DistributionMapping dm{ba};
    int num_grow = 0;
    amrex::MultiFab density(ba, dm, 1, num_grow);

    auto rho_arrs = density.arrays();
          amrex::ParallelFor(density, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {


  AMREX_ALWAYS_ASSERT(m2c::Length(1.0) == 100.0);
  AMREX_ALWAYS_ASSERT(c2m::Length(1.0, 2) == 1.0e-4);
  AMREX_ALWAYS_ASSERT(m2c::Mass(c2m::Mass(1.0)) == 1.0);
  AMREX_ALWAYS_ASSERT(m2c::Mu(1.0) == m2c::Mass(1.0) / m2c::Length(1.0));
  //amrex::Print() << "Unit Conversion Tests Passed" << std::endl;

  // Test locate function
  amrex::Array<amrex::Real, 5> sorted_data{0.0, 0.2, 0.5, 0.9, 1.0};
  int length = sorted_data.size();
  int idxlo = -1;
  amrex::Real x = -50.0;
  pele::physics::utilities::locate(sorted_data.data(), length, x, idxlo);
  AMREX_ALWAYS_ASSERT(idxlo == 0);
  x = 1.0;
  pele::physics::utilities::locate(sorted_data.data(), length, x, idxlo);
  AMREX_ALWAYS_ASSERT(idxlo == length - 1);
  x = 0.25;
  pele::physics::utilities::locate(sorted_data.data(), length, x, idxlo);
  AMREX_ALWAYS_ASSERT(idxlo == 1);
  //amrex::Print() << "Locate Tests Passed" << std::endl;

  // Test rectangle circle intersection area
  // Quarter circle intersection
  amrex::Real x0 = 0.0;
  amrex::Real x1 = 1.0;
  amrex::Real y0 = 0.0;
  amrex::Real y1 = 1.0;
  amrex::Real cx = 1.0;
  amrex::Real cy = 1.0;
  amrex::Real r = 1.0;
  amrex::Real A1 = pele::physics::utilities::rectangle_circle_intersection_area(
    x0, x1, y0, y1, cx, cy, r);
  AMREX_ALWAYS_ASSERT(std::abs(A1 - 0.25 * M_PI) < 1e-14);
  // No intersection
  cx = 3.0;
  cy = 3.0;
  amrex::Real A2 = pele::physics::utilities::rectangle_circle_intersection_area(
    x0, x1, y0, y1, cx, cy, r);
  AMREX_ALWAYS_ASSERT(std::abs(A2) < 1e-14);
  // Square covered by circle
  r = 20.0;
  amrex::Real A3 = pele::physics::utilities::rectangle_circle_intersection_area(
    x0, x1, y0, y1, cx, cy, r);
  AMREX_ALWAYS_ASSERT(std::abs(A3 - 1.0) < 1e-14);
  //amrex::Print() << "Intersection Area Tests Passed" << std::endl;

  rho_arrs[box_no]( i,j,k) = A1 + A2 + A3 + idxlo + m2c::Length(1.0);
        });
  amrex::Finalize();
  return 0;
}
