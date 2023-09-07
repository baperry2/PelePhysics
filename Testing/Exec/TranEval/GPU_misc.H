#ifndef GPU_MISC_H
#define GPU_MISC_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <cmath>

#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
initialize_data(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& spec,
  amrex::Array4<amrex::Real> const& temp,
  amrex::Array4<amrex::Real> const& rho,
  amrex::GeometryData const& geomdata) noexcept
{
  const amrex::Real* plo = geomdata.ProbLo();
  const amrex::Real* phi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();

  amrex::Real dTemp = 5.0;
  amrex::Real dRho = 0.005;
  amrex::Real y = plo[1] + (j + 0.5) * dx[1];
  amrex::Real x = plo[0] + (i + 0.5) * dx[0];
  amrex::Real pi = 3.1415926535897932;
  amrex::GpuArray<amrex::Real, 3> L;
  amrex::GpuArray<amrex::Real, 3> P;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_lo;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_hi;

  for (int n = 0; n < 3; n++) {
    L[n] = phi[n] - plo[n];
    P[n] = L[n] / 4.0;
  }
  for (int n = 0; n < NUM_SPECIES; n++) {
    Y_lo[n] = 0.0;
    Y_hi[n] = 1.0 / NUM_SPECIES;
  }
  Y_lo[0] = 1.0;

  // T, Yk, rho
#if (AMREX_SPACEDIM == 1)
  temp(i, j, k) = 500.0;
  rho(i, j, k) = 0.01;
#else
  temp(i, j, k) = 500.0 + dTemp * std::sin(2.0 * pi * (y - plo[1]) / P[1]);
  rho(i, j, k) = 0.01 + dRho * std::sin(2.0 * pi * (y - plo[1]) / P[1]);
#endif
  for (int n = 0; n < NUM_SPECIES; n++) {
    spec(i, j, k, n) = Y_lo[n] + (Y_hi[n] - Y_lo[n]) * (x - plo[0]) / L[0];
  }
  // corr Yk
  amrex::Real dummy = 0.0;
  for (int n = 0; n < NUM_SPECIES - 1; n++) {
    dummy = dummy + spec(i, j, k, n);
  }
  spec(i, j, k, NUM_SPECIES - 1) = 1.0 - dummy;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
initialize_data_linear(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& spec,
  amrex::Array4<amrex::Real> const& temp,
  amrex::Array4<amrex::Real> const& rho,
  amrex::GeometryData const& geomdata,
  amrex::Real min_T,
  amrex::Real max_T,
  amrex::Real min_P,
  amrex::Real max_P,
  amrex::Real min_Y_H2O,
  amrex::Real max_Y_H2O) noexcept
{
  const amrex::Real* plo = geomdata.ProbLo();
  const amrex::Real* phi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();

  amrex::IntVect idx = {AMREX_D_DECL(i,j,k)};
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> alpha;

  for (int n = 0; n < 3; n++) {
    alpha[n] = ((idx[n]+0.5)*dx[n])  / (phi[n] - plo[n]);
  }

#ifdef H2O_ID
  int idx_0 = H2O_ID;
#else
  int idx_0 = 0;
#endif
#ifdef CO2_ID
  int idx_1 = CO2_ID;
#else
  int idx_1 = NUM_SPECIES-1;
#endif
  if (idx_0 == idx_1) {
    amrex::Abort("initialize_data_linear: Can't set up compositions properly");
  }

  amrex::GpuArray<amrex::Real, NUM_SPECIES> mf_arr = {0.0};

  // T, Yk, rho
#if (AMREX_SPACEDIM == 1)
  temp(i, j, k) = min_T + (max_T - min_T) * alpha[0];
  amrex::Real press = min_P;
  mf_arr[idx_0] = min_Y_H2O;
  mf_arr[idx_1] = 1.0 - mf_arr[idx_0];
#elif (AMREX_SPACEDIM == 2)
  temp(i, j, k) = min_T + (max_T - min_T) * alpha[0];
  amrex::Real press = min_P + (max_P - min_P) * alpha[1];
  mf_arr[idx_0] = min_Y_H2O;
  mf_arr[idx_1] = 1.0 - mf_arr[idx_0];
#else
  temp(i, j, k) = min_T + (max_T - min_T) * alpha[0];
  amrex::Real press = min_P + (max_P - min_P) * alpha[1];
  mf_arr[idx_0] = min_Y_H2O + (max_Y_H2O - min_Y_H2O) * alpha[2];
  mf_arr[idx_1] = 1.0 - mf_arr[idx_0];
#endif
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2R(press, mf_arr.data(), temp(i, j, k), rho(i,j,k));
  for (int n = 0; n < NUM_SPECIES; ++n) {
    spec(i,j,k,n) = mf_arr[n];
  }
}
#endif
