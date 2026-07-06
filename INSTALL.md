# Documentation for Installation 

## Dependency & Environment

This FE library requires the following software and libraries:
- Fortran compiler supporting Fortran 2008
- MPI library
- LAPACK
- NetCDF
- SCALE library (https://scale.riken.jp)

For visualization, a Python environment with several Python libraries, such as Xarray and netCDF4, is required.
Note that Matplotlib version 3.8 or later is assumed.

We have confirmed that our codes can be built successfully in the following environments:
- Ubuntu 24.04 LTS (with GNU compiler)
  - GNU Fortran 13.2.0
  - Open MPI 4.1.6
  - LAPACK 3.12.0
  - NetCDF 4.9.2
  - SCALE library 5.5.5
- Ubuntu 22.04 LTS (with Intel oneAPI HPC Toolkit 2025.0)
  - Intel® Fortran Compiler
  - Intel MPI Library
  - Intel oneAPI Math Kernel Library
  - NetCDF 4.8.0
  - SCALE library develop version
- Ubuntu 22.04 LTS (with GNU compiler)
  - GNU Fortran 11.4.0
  - Open MPI 4.1.2
  - LAPACK 3.10.0
  - NetCDF 4.8.1
  - SCALE library 5.5.5
- Ubuntu 22.04 LTS (with NVIDIA HPC SDK 26.1 for OpenACC/GPU execution)
  - NVIDIA Fortran compiler (nvfortran)
  - NVIDIA Performance Libraries (NVPL)
  - CUDA-aware Open MPI, or another MPI library compatible with nvfortran
  - NetCDF 4.8.1
- macOS Tahoe
  - GNU Fortran 15.2.0
  - Open MPI 5.0.8
  - LAPACK 3.12.0
  - NetCDF 4.9.3
  - SCALE library develop version

Our codes have also been verified on Fugaku and Odyssey with Fujitsu compilers.

## Build FE-library

1. Preparation
  - Set the SCALE_FE_SYS environment variable. See the sysdep directory for available system-dependent settings.

    `% export SCALE_FE_SYS=MacOSX-gnu-ompi`   (for example)

  - Set the directory that contains SCALE library.

    `% export SCALE="~/workspace/scale-5.5.5/"`   (for example)

  - To enable thread parallelization with OpenMP, set the following variable:
    `% export SCALE_ENABLE_OPENMP=T`

  - Set the directory that contains the NetCDF library, if necessary.

    `% export NETCDF="/ap/netcdf4-fortran/4.7.3/"`   (for example)

  - If a development version of the SCALE library is used, set the following variable:

    `% export SCALE_DEVELOP=T`


2. Build the library in the FElib source directory.

    `% cd rootdir/FElib/src/`

    `% make`

## Compile and Run Simple Sample Programs

For example, in the case of sample/advect1d, 
 
  `% cd rootdir/sample/advect1d/`

  `% make`

  `% make run`


## Compile and Run Atmospheric Models

 To build the three-dimensional nonhydrostatic atmospheric model (SCALE-DG) and conduct an idealized test case, such as the density-current case, run the following commands:
 
  `% cd rootdir/model/atm_nonhydro3d/test/case/density_current`

  `% make`

  `% make run`

 The visualize directory includes several scripts using Python libraries.
 To visualize the simulation results using these scripts, run:

  `% make vis`

## Switches for Compilation

The following options are available as environment variables:

  `% export OPTION_NAME=T`

Switch options are activated only when the value is set to T.

- SCALE_ENABLE_OPENMP   : Enable thread parallelization with OpenMP
- SCALE_ENABLE_OPENACC  : Enable OpenACC directives for GPU execution in supported programs and components (experimental)
- SCALE_IGNORE_SRCDEP   : Ignore source dependency during compilation
 