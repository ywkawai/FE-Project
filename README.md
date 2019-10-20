FE-project 

=======================================================================================

Target
-----------------------------------------------------------------------------------------
FE project provides a library for some discontinuous Galerkin methods (DGMs), 
and some sample programs with the library for DGMs. 

Sample programs
-----------------------------------------------------------------------------------------
* 1-dimensional linear advection problems
* 2-dimensional linear advection problems in a rectangle domain
* a 2-dimensional nonhydrostatic atmospheric model

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

3. compile and run 

 For simple sample programs such as sample/advect1d or sample/advect2d, 
 e.g.,
 
 `% cd rootdir/sample/advect1d/`

 `% make`

 If you want to build a 2-dimensional nonhydrostatic atmospheric model, 
 and conduct a idealized test case, such as density current, using it, 
 
 `% cd rootdir/sample/atm_nonhydro2d/test/case/density_current`

 `% make`
 `% make run`

 In the directory of 'visualize', some python scripts with matplotlib 
 are prepared for visualizing simulation results. 

 Gallery of results
 ----------------------------------------------------------------------------------------
 * a 2-dimensional nonhydrostatic atmospheric model
   ![A simulation of density current](https://raw.github.com/wiki/ywkawai/FE-projet/gallery/atm_nohydro2d/density_current/density_current.gif)

 For more animations of simulation results, please see the Wiki of FE-project on github or 'FE-project gallery' channel on YouTube (url: https://www.youtube.com/channel/UCO17OQtKHwkkQwmHD9y9mQg/featured). 

