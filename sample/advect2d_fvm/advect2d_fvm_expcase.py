#!/usr/bin/env python
#------------------------------------------------------
import textwrap
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../auxiliary'))
from experiment import Experiment, ExpCase

class Advect2DCase(ExpCase):
  def __init__(self, exp, path_from_exptopdir, name, neX, neY, flux_scheme, tint_scheme, duration_sec, dt_sec):
    self._neX = neX; self._neY = neY
    self._flux_scheme = flux_scheme
    self._tint_scheme = tint_scheme
    super(Advect2DCase, self).__init__( exp, path_from_exptopdir, name, "test.conf",   
                                        duration_sec, dt_sec, 
                                       "test_advect2d_fvm", 1, 1)
  
  def set_initShape(self, iniShapeName, iniShapeParams):
    self._iniShapeName = iniShapeName
    self._iniShapeParams = iniShapeParams

  def set_velType(self, velTypeName, velTypeParams):
    self._velTypeName = velTypeName
    self._velTypeParams = velTypeParams

  def generate_configfile(self): 

    ishape_params = "D0, ".join(map(str, self._iniShapeParams)) + "D0"
    veltype_params = "D0, ".join(map(str, self._velTypeParams)) + "D0"
    conf = []
    nstep_eval_error = int(0.1/self._dt_sec)
    conf.append( textwrap.dedent('''\
        &PARAM_TEST
          NeGX             = {neX}, 
          GXHALO           = 2, 
          NeGY             = {neY}, 
          GYHALO           = 2, 
          InitShapeName    = '{iniShapeName}',
          InitShapeParams  = {ishape_params},
          VelTypeName    = '{velTypeName}',
          VelTypeParams  = {veltype_params},        
          nstep_eval_error = {nstep_eval_error}, 
          FLUX_SCHEME_TYPE = {flux_scheme_type}, 
          TINTEG_SCHEME_TYPE = 'RK_TVD_3',  
        /
    ''').format( \
      neX=self._neX, neY=self._neY, flux_scheme_type=self._flux_scheme, 
      iniShapeName=self._iniShapeName, ishape_params=ishape_params, 
      velTypeName=self._velTypeName, veltype_params=veltype_params, nstep_eval_error=nstep_eval_error) )
    
    conf.append( textwrap.dedent('''\
        &PARAM_TIME
          TIME_DURATION                 = {duration}D0, 
          TIME_DURATION_UNIT            = 'SEC', 
          TIME_DT                       = {dt}D0, 
          TIME_DT_UNIT                  = 'SEC', 
        /
    ''').format(duration=self._duration_sec, dt=self._dt_sec ) )
  
    conf.append( textwrap.dedent('''\
        &PARAM_FILE_HISTORY
        FILE_HISTORY_DEFAULT_BASENAME  = "history",
        FILE_HISTORY_DEFAULT_TINTERVAL = 0.1D0,
        FILE_HISTORY_DEFAULT_TUNIT     = "SEC",
        FILE_HISTORY_DEFAULT_TAVERAGE  = .false.,
        FILE_HISTORY_DEFAULT_DATATYPE  = "REAL4",
        FILE_HISTORY_OUTPUT_STEP0      = .true.,
        /
        &HISTORY_ITEM name='q'         /
        &HISTORY_ITEM name='qexact'    /      
    ''') )
    with open(self._path+"/test.conf", mode='w') as f:
      f.write( "\n".join(conf) )


def exp_prepair(TOP_DIR, Ne_List, FluxScheme_List, DELT_List, ORG_PE_PATH, expcase_params_dict): 
  exp = Experiment(TOP_DIR, "advect2d_fvm", ORG_PE_PATH)

  for fs_i, flux_scheme in enumerate(FluxScheme_List):
    for ne_i, ne in enumerate(Ne_List):
      for expcase_prefix, params in expcase_params_dict.items():
        iniShapeName = params[0]; iniShapeParams = params[1] 
        velTypeName = params[2]; velTypeParams = params[3]

        path = "I"+iniShapeName+"/V"+velTypeName+"/"+"Ne"+str(ne)+str(flux_scheme)
        case_name = expcase_prefix + "_Ne"+str(ne)+str(flux_scheme)

        case = Advect2DCase(exp, path, case_name, ne, ne, flux_scheme, 'RK_TVD_3', 5.0, DELT_List[fs_i][ne_i])
        case.set_initShape(iniShapeName, iniShapeParams)
        case.set_velType(velTypeName, velTypeParams)
        exp.add_expcase(case)

  exp.prepair()

