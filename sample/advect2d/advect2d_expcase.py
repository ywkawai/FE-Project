#!/usr/bin/env python
#------------------------------------------------------
import textwrap
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../auxiliary'))
from experiment import Experiment, ExpCase

class Advect2DCase(ExpCase):
  def __init__(self, exp, path_from_exptopdir, name, neX, neY, polyOrder, tint_scheme, duration_sec, dt_sec):
    self._neX = neX; self._neY = neY
    self._polyOrder = polyOrder
    self._tint_scheme = tint_scheme
    super(Advect2DCase, self).__init__( exp, path_from_exptopdir, name, "test.conf",   
                                        duration_sec, dt_sec, 
                                       "test_advect2d", 1, 1)
  
  def set_initShape(self, iniShapeName, iniShapeParams, iniGalerkinProjFlag):
    self._iniShapeName = iniShapeName
    self._iniShapeParams = iniShapeParams
    self._iniGalerkinProjFlag = iniGalerkinProjFlag

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
          NeGY             = {neY}, 
          PolyOrder        = {polyOrder},
          InitShapeName    = '{iniShapeName}',
          InitShapeParams  = {ishape_params},
          InitCond_GalerkinProjFlag = {iniGalerkinProjFlag}, 
          VelTypeName    = '{velTypeName}',
          VelTypeParams  = {veltype_params},        
          nstep_eval_error = {nstep_eval_error}, 
          TINTEG_SCHEME_TYPE = 'RK_TVD_3',  
        /
    ''').format( \
      neX=self._neX, neY=self._neY, polyOrder=self._polyOrder, 
      iniShapeName=self._iniShapeName, ishape_params=ishape_params, iniGalerkinProjFlag=self._iniGalerkinProjFlag, 
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


def exp_prepair(TOP_DIR, Ne_List, PolyOrder_List, DELT_List, ORG_PE_PATH, expcase_params_dict, iniGalerkinProjFlag): 
  exp = Experiment(TOP_DIR, "advect2d", ORG_PE_PATH)

  for p_i, p in enumerate(PolyOrder_List):
    for ne_i, ne in enumerate(Ne_List):
      print(p_i, ne_i)
      for expcase_prefix, params in expcase_params_dict.items():
        iniShapeName = params[0]; iniShapeParams = params[1] 
        velTypeName = params[2]; velTypeParams = params[3]

        path = "I"+iniShapeName+"/V"+velTypeName+"/"+"Ne"+str(ne)+"P"+str(p)
        case_name = expcase_prefix + "_Ne"+str(ne)+"P"+str(p)

        case = Advect2DCase(exp, path, case_name, ne, ne, p, 'RK_TVD_3', 5.0, DELT_List[p_i][ne_i])
        case.set_initShape(iniShapeName, iniShapeParams, iniGalerkinProjFlag)
        case.set_velType(velTypeName, velTypeParams)
        exp.add_expcase(case)

  exp.prepair()

