# AMReX
DIM        = 3
PRECISION  = DOUBLE
PROFILE    = FALSE
VERBOSE    = FALSE
DEBUG      = FALSE

# Compiler
COMP       = gnu
FCOMP      = gfortran
USE_MPI    = TRUE
USE_OMP    = FALSE
USE_CUDA   = FALSE
USE_HIP    = FALSE
USE_LIBTORCH = FALSE

# PelePhysics
FUEGO_GAS  = TRUE
USE_MANIFOLD = TRUE
TINY_PROFILE = TRUE

# define the location of the PELE_PHYSICS top directory
PELE_PHYSICS_HOME    ?= ../../..
# TODO FIXME DELETE HACK HACK
AMREX_HOME = /Users/bperry/Software/amrex
SUNDIALS_HOME = /Users/bperry/Software/sundials

# this flag activates the subcycling mode in the D/Cvode routines
DEFINES  += -DMOD_REACTOR

# Activates use of SUNDIALS
ifeq ($(USE_CUDA), TRUE)
  PELE_USE_KLU = FALSE
else
  ifeq ($(USE_HIP), TRUE)
    PELE_USE_KLU = FALSE
  else
    PELE_USE_KLU = FALSE
  endif
endif

ifeq ($(FUEGO_GAS), TRUE)
  Eos_Model       = Fuego
  Chemistry_Model = grimech30-noArN
  Transport_Model = Constant
else ifeq ($(USE_MANIFOLD), TRUE)
  Eos_Model       = Manifold
  Chemistry_Model = Null
  Transport_Model = Manifold
  Manifold_Dim    = 2
else
  Eos_Model       = GammaLaw
  Chemistry_Model = Null
  Transport_Model = Constant
endif

Bpack   := ./Make.package
Blocs   := .

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
