#!/bin/bash
#----------- install module -------------
module load netcdf-c
module load netcdf-fortran
module load parallel-netcdf
module load phdf5

module load python/3.8.12
module load py-numpy
module load py-matplotlib

#---- Echo my work dir ----
echo "My work directory: ${WORK_DIR}"

#------------ setup environment variables ---------
export SCALE=${WORK_DIR}/scale-5.4.5
export SCALE_FE=${WORK_DIR}/FE-Project
export SCALE_ENABLE_OPENMP=T
export SCALE_FE_SYS=FLOW
export SCALE_SYS=FLOW
export FE_SRC=$SCALE_FE/FElib/src

#---------- print out enviroment variables -----
echo $SCALE
echo $SCALE_FE
echo $SCALE_FE_SYS
echo $SCALE_SYS

