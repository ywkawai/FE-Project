FE-project 

=======================================================================================

Target
-----------------------------------------------------------------------------------------
FE project provides a library for some discontinuous Galerkin methods (DGMs), 
and some sample programs with the library for DGMs. 

Sample programs
-----------------------------------------------------------------------------------------
* 1-dimensional linear advection problem
* 2-dimensional linear advection problem in a rectangle domain

Dependency
----------------------------------------------------------------------------------------
This FE library requires SCALE library (http://r-ccs-climate.riken.jp/scale/ja/index.html). 

Build
----------------------------------------------------------------------------------------
1. preparation
  - set SCALE_FE_SYS environmental variable (see the sysdef directory)

  `% export SCALE_FE_SYS=MacOSX-gnu-ompi`   (for example)

  - set a directory in which SCALE library is contained

  `% export SCALE="~/workspace/scale-5.3.3/"`   (for example)

  - set a directory in which NetCDF library is contained (if necessary).
    
  `% export NETCDF="/ap/netcdf4-fortran/4.4.2/"`   (for example)

2. build the library in the directory of FElib

 `% cd rootdir/FElib/src/`

 `% make`

3. compile and run a sample program in the directory of sample
 e.g.,
 
 `% cd rootdir/sample/advect1d/`

 `% make`
 