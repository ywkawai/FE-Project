################################################################################
#
# for Fugaku
#
################################################################################
#PJM --name "GW"
#PJM -L rscgrp=regular-o
#PJM -L node=6
#PJM -g ga41
#PJM --rsc-list "elapse=00:05:00"
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

SCALE_DG_INIT_BIN=../scale-dg_init
SCALE_DG_BIN=../scale-dg
SCALE_DG_REGRID_BIN=../../../../../../bin/regrid_tool
  
#llio_transfer ${SCALE_DG_INIT_BIN} ${SCALE_DG_BIN} *.conf

mpiexec -np 24 -of-proc ./output.%j/%/1000r/stdout -oferr-proc ./output.%j/%/1000r/stderr \
  ${SCALE_DG_INIT_BIN} init.conf || exit 1    
    
mpiexec -np 24 -ofout-proc ./output.%j/%/1000r/stdout -oferr-proc ./output.%j/%/1000r/stderr \
  ${SCALE_DG_BIN} run.conf || exit 1  

#llio_transfer --purge ${SCALE_DG_INIT_BIN} ${SCALE_DG_BIN} *.conf    
  