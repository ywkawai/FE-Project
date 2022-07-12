# Documentation for installation 

## Dependency & Environment

This FE library requires following libraries: 
  - Fortran compiler supporting Fortran 2008
  - MPI library
  - LAPACK
  - NetCDF
  - SCALE library (https://scale.riken.jp). 

We confirm that building our codes has succeeded in the following environments:
  - Ubuntu 20.04 LTS (for case of GNU compiler)
    - GNU Fortran 9.3.0
    - OpenMPI 4.0.3
    - LAPACK 3.9.0
    - NetCDF 4.7.3
    - SCALE library 5.4.5
  - Ubuntu 20.04 LTS (for case of Intel oneAPI HPC Toolkit 2022.1)
    - Intel® Fortran Compiler
    - Intel MPI Library
    - Intel oneAPI Math Kernel Library
    - NetCDF 4.8.0
    - SCALE library develop version
  - macOS Monterey
    - GNU Fortran 11.3.0
    - OpenMPI 4.1.3
    - LAPACK 3.10.1
    - NetCDF 4.8.1
    - SCALE library develop version

Our codes are also verified in Oakbridge-CX (Intel compiler) and Fugaku (Fujitsu compiler). 

## Build FE-library

1. preparation
  - set SCALE_FE_SYS environmental variable (see the sysdef directory)

  `% export SCALE_FE_SYS=MacOSX-gnu-ompi`   (for example)

  - set a directory in which SCALE library is contained

  `% export SCALE="~/workspace/scale-5.4.5/"`   (for example)

  - If you use the develop version of SCALE library, set a variable as

  `% export SCALE_DEVELOP=T`

  - set a directory in which NetCDF library is contained (if necessary).

  `% export NETCDF="/ap/netcdf4-fortran/4.7.3/"`   (for example)

2. build the FE-library

 `% cd rootdir/`

 `% make`

## Compile and run simple sample programs

 In rootdir/sample/ directory, there are simple sample programs. To compile and run it, for example, in the case of sample/advect1d, 
 
 `% cd rootdir/sample/advect1d/`

 `% make`

 `% make run`


## Build atmospheric models and Perform numerical experiments

 If you would like to conduct idealized test cases, such as density current using a three-dimensional nonhydrostatic model, 
 
 `% cd rootdir/model/atm_nonhydro3d/test/case/density_current`

 `% make`

 If this procedure succeeded, three binary files (scale-dg, scale-dg_init, scale-dg_pp) are generated. 

  By executing 

 `% make jobshell`

  we can prepare a job script named run.sh.

  Finally, 

 `% sh run.sh`

 In the directory of 'visualize', some python scripts with matplotlib are prepared for visualizing simulation results. 
