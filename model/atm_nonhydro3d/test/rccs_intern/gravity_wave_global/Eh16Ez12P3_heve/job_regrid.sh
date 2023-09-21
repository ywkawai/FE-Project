################################################################################
#
# for Fugaku
#
################################################################################
#PJM --name "REGRID"
#PJM -L rscgrp=regular-o
#PJM -L node=16
#PJM -g gz06
#PJM --rsc-list "elapse=00:30:00"
#PJM --mpi "max-proc-per-node=4"
#PJM -S

module load netcdf
module load netcdf-fortran
module load hdf5
module load pnetcdf

#export XOS_MMM_L_ARENA_FREE=1
export FORT90L="-Wl,-T"
#export OMPI_MCA_plm_ple_memory_allocation_policy=bind_local
export PLE_MPI_STD_EMPTYFILE="off"
export PARALLEL=12
export OMP_NUM_THREADS=12
#export fu11bf=1

SCALE_DG_REGRID_BIN=../../../../../../bin/regrid_tool
  
#llio_transfer ${SCALE_DG_REGRID_BIN} *.conf

mkdir -p outdata
mpiexec -np 64 -of-proc ./output.%j/%/1000r/stdout -oferr-proc ./output.%j/%/1000r/stderr \
  ${SCALE_DG_REGRID_BIN} regrid.conf || exit 1    
    
#llio_transfer --purge ${SCALE_DG_REGRID_BIN} *.conf    
  