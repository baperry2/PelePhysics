#!/usr/bin/env bash

CHEMINP=CH4_ion_LuSk.dat
THERMINP=therm.dat
FINALFILE=mechanism.cpp

source ../mechanism_builder.sh

${FUEGO_PYTHON} "${FMC}" -mechanism=${CHEMINP} -thermo=${THERMINP} -name=${FINALFILE}

echo "Compiling ${FINALFILE}..."
