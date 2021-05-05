# Documentation for installation 

## Dependency & Environment

This FE library requires following libraries: 
  - Fortran compiler supporting Fortran 2008. 
  - MPI library
  - LAPACK
  - NetCDF
  - SCALE library (https://scale.riken.jp). 

We confirm that building our codes has succeeded in the following environments:
  - macOS Mojave
    - GNU Fortran 10.2.0
    - OpenMPI 4.1.1
    - LAPACK 3.9.1
    - NetCDF 4.7.4
    - SCALE library 5.4.3
  - Ubuntu 20.04 LTS
    - GNU Fortran 9.3.0
    - OpenMPI 4.0.3
    - LAPACK 3.7.1
    - NetCDF 4.7.3
    - SCALE library 5.4.3

Our codes are verified in Oakbridge-CX (intel compiler) and Fugaku. 

## Build FE-library 

1. preparation
  - set SCALE_FE_SYS environmental variable (see the sysdef directory)

  `% export SCALE_FE_SYS=MacOSX-gnu-ompi`   (for example)

  - set a directory in which SCALE library is contained

  `% export SCALE="~/workspace/scale-5.4.3/"`   (for example)

  - If you use the develop version of SCALE library, set a variable as

  `% export SCALE_DEVELOP=T`

  - set a directory in which NetCDF library is contained (if necessary).

  `% export NETCDF="/ap/netcdf4-fortran/4.7.3/"`   (for example)

2. build the library in the directory of FElib

 `% cd rootdir/FElib/src/`

 `% make`

## Compile and run simple sample programs

 For example, in the case of sample/advect1d, 
 
 `% cd rootdir/sample/advect1d/`

 `% make`

## Compile and run dynamical core

 If you want to build a 3-dimensional nonhydrostatic atmospheric model, 
 and conduct a idealized test case, such as density current, using it, 
 
 `% cd rootdir/model/atm_nonhydro3d/test/case/density_current`

 `% make`

 `% make run`

 In the directory of 'visualize', some python scripts with matplotlib 
 are prepared for visualizing simulation results. 
