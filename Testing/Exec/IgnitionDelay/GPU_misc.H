#ifndef GPU_MISC_H
#define GPU_MISC_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <cmath>

#include <PelePhysics.H>

#ifdef USE_MANIFOLD_EOS
typedef pele::physics::ManFuncParams::ManFuncData ManFuncData;
typedef pele::physics::NNFuncParams::NNFuncData NNFuncData;
#endif

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
initialize_data(
  int i,
  int j,
  int k,
  int fuel_id,
  int o2_id,
  int n2_id,
  int energy_type,
  amrex::Array4<amrex::Real> const& rhoY,
  amrex::Array4<amrex::Real> const& frcExt,
  amrex::Array4<amrex::Real> const& rhoE,
  amrex::Array4<amrex::Real> const& frcEExt,
  amrex::Real t0,
  amrex::Real equiv_ratio,
  amrex::Real p,
  amrex::GeometryData const& geomdata,
  const pele::physics::eos::EosParm<
  pele::physics::PhysicsType::eos_type>* eosparm) noexcept
{
  const amrex::Real* plo = geomdata.ProbLo();
  const amrex::Real* phi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();

  amrex::Real Temp_lo = 2000.0;
  amrex::Real Temp_hi = 2500.0;
  amrex::Real dTemp = 100.0;
  amrex::Real pressure = p;
  amrex::Real density = 0.0;
  amrex::Real energy = 0.0;
  amrex::Real temp = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> X;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y;

  amrex::Real y = plo[1] + (j + 0.5) * dx[1];
  // amrex::Real x = plo[0] + (i+0.5)*dx[0];
  amrex::Real pi = 3.1415926535897932;
  amrex::GpuArray<amrex::Real, 3> L;
  amrex::GpuArray<amrex::Real, 3> P;
  
  auto eos = pele::physics::PhysicsType::eos(eosparm);

  for (int n = 0; n < 3; n++) {
    L[n] = phi[n] - plo[n];
    P[n] = L[n] / 4.0;
  }
  // Y
  for (int n = 0; n < NUM_SPECIES; n++) {
    X[n] = Y[n] = 0.0;
  }
  // X[O2_ID]   = 0.7;
  // X[fuel_id] = 0.3;
#ifdef USE_MANIFOLD_EOS

  // NOTE: Assumes fuel is methane!
  amrex::Real xo2, yo2, xn2, yn2, xfuel, yfuel, xwinv;
  xo2 = 1.0 / (equiv_ratio / 2.0 + 1 + 0.79 / 0.21);
  xn2 = 0.79 * xo2 / 0.21;
  xfuel = 1.0 - xo2 - xn2;
  xwinv = 1.0 / (xo2*31.998000 + xn2*28.014000 + xfuel*16.043000);
  yo2 = xo2 * 31.998000 * xwinv;
  yn2 = xn2 * 28.014000 * xwinv;
  yfuel = xfuel * 16.043000 * xwinv;
  
  // Not clear what definition of molar fraction is for manifold variable
  
  if(eosparm->manf_data->manmodel == pele::physics::ManifoldModel::TABLE)
  {
    // Hardcode manifold variable definitions for table
    Y[0] = 0.0; // Progress variable
    Y[1] = yfuel; // Fuel fraction
  }
  else
  {
    NNFuncData* nnf_data = static_cast<NNFuncData*>(eosparm->manf_data);
    for(int i = 0; i < NUM_SPECIES-1; i++)
    {
      const int defn_start = i*nnf_data->Ncomb;
      
      Y[i] = nnf_data->comb_coeff[defn_start + o2_id]*yo2
           + nnf_data->comb_coeff[defn_start + n2_id]*yn2
           + nnf_data->comb_coeff[defn_start + fuel_id]*yfuel
           + nnf_data->man_bias[i];
    }
  }
  Y[NUM_SPECIES-1] = 1.0;
#else
  X[O2_ID] = 1.0 / (equiv_ratio / 2.0 + 1 + 0.79 / 0.21);
  X[N2_ID] = 0.79 * X[O2_ID] / 0.21;
  X[fuel_id] = 1.0 - X[O2_ID] - X[N2_ID];
  eos.X2Y(&X[0], &Y[0]);
#endif

  // T
  temp = t0;
  // temp =  Temp_lo + (Temp_hi-Temp_lo)*y/L[1] + dTemp *
  // std::sin(2.0*pi*y/P[1]);
#ifdef USE_MANIFOLD_EOS
  eos.PYT2R(pressure, &Y[0], temp, density);
#else
  if (energy_type == 1) {
    // get rho and E
    eos.PYT2RE(pressure, &Y[0], temp, density, energy);
  } else {
    // get rho and H
    eos.PYT2R(pressure, &Y[0], temp, density);
    eos.TY2H(temp, &Y[0], energy);
  }
#endif
  // Fill vect
  for (int n = 0; n < NUM_SPECIES; n++) {
    rhoY(i, j, k, n) = Y[n] * density;
    frcExt(i, j, k, n) = 0.0;
  }
  rhoY(i, j, k, NUM_SPECIES) = temp;

  rhoE(i, j, k) = energy * density;

  frcEExt(i, j, k) = 0.0;
}

#endif
