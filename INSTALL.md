# Documentation for installation 

## Dependency & Environment

This FE library requires following libraries: 
  - Fortran compiler supporting Fortran 2008
  - MPI library
  - LAPACK
  - NetCDF
  - SCALE library (https://scale.riken.jp). 

We confirm that building our codes has succeeded in the following environments:
  - Ubuntu 24.04 LTS (for case of GNU compiler)
    - GNU Fortran 13.2.0
    - OpenMPI 4.1.6
    - LAPACK 3.12.0
    - NetCDF 4.9.2
    - SCALE library 5.5.4
  - Ubuntu 22.04 LTS (for case of Intel oneAPI HPC Toolkit 2025.0)
    - Intel® Fortran Compiler
    - Intel MPI Library
    - Intel oneAPI Math Kernel Library
    - NetCDF 4.8.0
    - SCALE library develop version
  - macOS Sequoia
    - GNU Fortran 13.2.0
    - OpenMPI 5.0.3
    - LAPACK 3.12.0
    - NetCDF 4.9.2
    - SCALE library develop version

Our codes are also verified in Fugaku and Odyssey (Fujitsu compiler). 

## Build FE-library

1. Preparation
  - Set SCALE_FE_SYS environmental variable (see the sysdef directory)

  `% export SCALE_FE_SYS=MacOSX-gnu-ompi`   (for example)

  - Set a directory in which SCALE library is contained

  `% export SCALE="~/workspace/scale-5.5.4/"`   (for example)

  - If you would like to enable a thread parallelization with OpenMP, set a variable as 

  `% export SCALE_ENABLE_OPENMP=T`

  - Set a directory in which a NetCDF library is contained (if necessary).

  `% export NETCDF="/ap/netcdf4-fortran/4.7.3/"`   (for example)

  - If need, to indicate that a developing version of SCALE library is used, set a variable as

  `% export SCALE_DEVELOP=T`


2. Build the library in the directory of FElib

 `% cd rootdir/FElib/src/`

 `% make`

## Compile and run simple sample programs

 For example, in the case of sample/advect1d, 
 
 `% cd rootdir/sample/advect1d/`

 `% make`

 `% make run`


## Compile and run atmospheric models

 If you want to build a three-dimensional nonhydrostatic atmospheric model, 
 and conduct an idealized test case, such as density current, using it, 
 
 `% cd rootdir/model/atm_nonhydro3d/test/case/density_current`

 `% make`

 `% make run`

 The directory of 'visualize' includes several scripts using Python libraries (e.g., xarray and matplotlib). 
 To visualize the simulation result using the scripts, 

  `% make vis`
