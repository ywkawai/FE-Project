import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../common'))
import batch_job_common


ptracer_shape="GAUSSIAN"
lonc=4.71238898038469 # 3PI/2
latc=0 
rh=1274244 # RPlanet / 5.0
initgp_porder=11
tau=1036800
  
SCALE_DG_BIN_PATH="../../.."
SCALE_DG_REGRID_BIN_PATH="../../../../../../../../bin"

#----------------------

def mkconf_init( conf_path,
                nprc, neh, nez, porder, 
                ptracer_shape, lonc, latc, rh, initgp_porder, 
                tau, alph): 
    conf_init_s = f"""#--- Configuration file for a test case of two-dimensional linear advection  -------
&PARAM_IO
 IO_LOG_BASENAME = 'init_LOG',
/
&PARAM_MKINIT
  initname = 'tracer_advection_global', 
/
&PARAM_RESTART
  OUTPUT_FLAG = .true., 
  OUT_BASENAME = 'init'
/
&PARAM_TIME
  TIME_STARTDATE       = 0000, 1, 1, 0, 0, 0,
  TIME_STARTMS         = 0.D0,
/
&PARAM_CONST
 CONST_OHM = 0D0, 
/
&PARAM_EXP
  FLOW_TYPE = 'SOLID_BODY_ROTATION_FLOW', 
  SOLID_BODY_ROT_TAU  = {tau}D0, 
  SOLID_BODY_ROT_ALPH = {alph}D0, 
  INIT_TRACER_PROF    = '{ptracer_shape}', 
  LONC                = {lonc}D0,   
  LATC                = {latc}D0, 
  {ptracer_shape}_rh         = {rh}D0, 
  InitGP_PolyOrder_h  = {initgp_porder}, 
  InitGP_PolyOrder_v  = {initgp_porder},
/
#** ATMOS ******************************************************
&PARAM_ATMOS
  ATMOS_MESH_TYPE = 'GLOBAL', 
  ACTIVATE_FLAG = .true., 
  ATMOS_DYN_DO  = .true.
/
&PARAM_ATMOS_MESH
  NLocalMeshPerPrc = 1, 
  Nprc             = {nprc}, 
  NeGX             = {neh},
  NeGY             = {neh},
  NeZ              = {nez},
  dom_zmin         = 0.0D0, 
  dom_zmax         = 10.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},
!  LumpedMassMatFlag = .true.,   
/
#** ATMOS / DYN ******************************************************
&PARAM_ATMOS_DYN
  EQS_TYPE         = "NONHYDRO3D_HEVE", 
  TINTEG_TYPE      = 'ERK_SSP_3s3o',    
/    
    """
    
    with open(conf_path, 'w') as f:
        f.write(conf_init_s)

#----------------

def mkconf_run( conf_path,
                nprc, neh, nez, porder, 
                ptracer_shape, lonc, latc, rh, initgp_porder, 
                tau, alph, 
                dt ): 
    conf_run_s = f"""#--- Configuration file for a test case of two-dimensional linear advection  -------
&PARAM_RESTART
  IN_BASENAME = "init_00000101-000000.000",
  OUTPUT_FLAG = .true., 
  OUT_BASENAME = 'restart'
/
&PARAM_TIME
  TIME_STARTDATE       = 0000, 1, 1, 0, 0, 0,
  TIME_STARTMS         = 0.D0,
  TIME_DURATION        = 12.0D0, 
  TIME_DURATION_UNIT   = 'DAY', 
  TIME_DT              = 1200.0D0, 
  TIME_DT_UNIT         = 'SEC', 
/
&PARAM_CONST
 CONST_OHM = 0D0, 
/
&PARAM_USER
  USER_do = .true., 
  polyOrderErrorCheck = {initgp_porder},
  LOG_STEP_INTERVAL   = 1,    
/
&PARAM_EXP
  FLOW_TYPE = 'SOLID_BODY_ROTATION_FLOW', 
  SOLID_BODY_ROT_TAU  = {tau}D0, 
  SOLID_BODY_ROT_ALPH = {alph}D0, 
  INIT_TRACER_PROF    = '{ptracer_shape}', 
  LONC                = {lonc}D0,   
  LATC                = {latc}D0, 
  {ptracer_shape}_rh         = {rh}D0, 
  InitGP_PolyOrder_h  = {initgp_porder}, 
  InitGP_PolyOrder_v  = {initgp_porder},
/
#** ATMOS ******************************************************
&PARAM_ATMOS
  ACTIVATE_FLAG       = .true., 
  ATMOS_MESH_TYPE     = 'GLOBAL',   
  TIME_DT             = 1200.0D0, 
  TIME_DT_UNIT        = 'SEC', 
  ATMOS_DYN_DO        = .true.
/
&PARAM_ATMOS_MESH
  NLocalMeshPerPrc = 1, 
  Nprc             = {nprc}, 
  NeGX             = {neh},
  NeGY             = {neh},
  NeZ              = {nez},
  dom_zmin         = 0.0D0, 
  dom_zmax         = 10.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},
!  LumpedMassMatFlag = .true.,   
/
&PARAM_ATMOS_VARS
  CHECK_RANGE = .true. ,
  CHECK_TOTAL = .false.,
/
#** ATMOS / DYN ******************************************************
&PARAM_ATMOS_DYN
  EQS_TYPE            = "NONHYDRO3D_HEVE", 
  ONLY_TRACERADV_FLAG       = .true., 
  TRACERADV_DISABLE_LIMITER = .true., 
  !-
  TINTEG_TYPE         = 'ERK_SSP_10s4o_2N', ! Dummy because ONLY_TRACERADV_FLAG = .true.
  TINTEG_TYPE_TRACER  = 'ERK_SSP_10s4o_2N', !  
  TIME_DT             = {dt}D0, 
  TIME_DT_UNIT        = 'SEC', 
  !-
  MODALFILTER_FLAG  = .false.,
  NUMDIFF_FLAG      = .false.,
/
&PARAM_ATMOS_DYN_BND
  btm_vel_bc   = 'SLIP', 
  top_vel_bc   = 'SLIP', 
  north_vel_bc = 'PERIODIC',
  south_vel_bc = 'PERIODIC', 
  east_vel_bc  = 'PERIODIC',
  west_vel_bc  = 'PERIODIC', 
  btm_thermal_bc   = 'ADIABATIC', 
  top_thermal_bc   = 'ADIABATIC', 
  north_thermal_bc = 'PERIODIC',
  south_thermal_bc = 'PERIODIC', 
  west_thermal_bc  = 'PERIODIC',
  east_thermal_bc  = 'PERIODIC',
/

#*** OUTPUT *******************************************
&PARAM_FILE_HISTORY
 FILE_HISTORY_DEFAULT_BASENAME  = "history",
 FILE_HISTORY_DEFAULT_TINTERVAL = 0.5D0, 
 FILE_HISTORY_DEFAULT_TUNIT     = "DAY",
 FILE_HISTORY_DEFAULT_TAVERAGE  = .false.,
 FILE_HISTORY_DEFAULT_DATATYPE  = "REAL4",
 FILE_HISTORY_OUTPUT_STEP0      = .true.,
/
&HISTORY_ITEM name='U'        /
&HISTORY_ITEM name='V'        /
&HISTORY_ITEM name='W'        /
&HISTORY_ITEM name='PTracer'  /

#*** Statistics *******************************************

&PARAM_MESHFIELD_STATISTICS
 use_globalcomm = .true.,
/
&PARAM_MONITOR 
  MONITOR_STEP_INTERVAL = 1 
/
&MONITOR_ITEM name='PTracer' /
    """
    
    with open(conf_path, 'w') as f:
        f.write(conf_run_s)

def mkconf_regrid( conf_path,
                nprc, neh, nez, porder, 
                regrid_nprcx, regrid_nprcy, 
                regrid_nex, regrid_ney, regrid_nez, 
                regrid_porder ): 
    conf_run_s = f"""#--- Configuration file for a test case of two-dimensional linear advection -------
&PARAM_IO
 IO_LOG_BASENAME = "regrid_LOG"
! IO_LOG_ALLNODE  = .true., 
/
&PARAM_REGRID_MESH
 in_MeshType  = "CUBEDSPHERE3D", 
 out_MeshType = "LONLAT3D", 
/
&PARAM_REGRID_INTERP_FIELD
  !- input --------------------
  in_basename="history",      
  vars = "U", "V", "PTracer",
  !out_tinterval = 5,
/
&PARAM_REGRID_FILE
  !-- output ----------------
  out_basename="./outdata/history", 
  out_UniformGrid=.false., 
/
&PARAM_REGRID_OPERATE_FIELD
  uvmet_conversion_flag = .true., 
/
&PARAM_REGRID_INMESH3D_CUBEDSPHERE
  NLocalMeshPerPrc = 1, 
  Nprc             = {nprc}, 
  NeGX             = {neh},
  NeGY             = {neh},
  NeGZ             = {nez},
  dom_zmin         = 0.0D0, 
  dom_zmax         = 10.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},
/
&PARAM_REGRID_OUTMESH3D_STRUCTURED
  NprcX       = {regrid_nprcx},       
  NeX         = {regrid_nex},           
  NprcY       = {regrid_nprcy}, 
  NeY         = {regrid_ney},    
  NeGZ        = {regrid_nez}, 
  PolyOrder_h = {regrid_porder}, 
  PolyOrder_v = {porder}, 
  dom_xmin    =   0.0D0, 
  dom_xmax    = 360.0D0,   
  dom_ymin    = -90.0D0, 
  dom_ymax    =  90.0D0, 
  dom_zmin    = 0.0D0, 
  dom_zmax    = 10.0D3,   
/
    """
    
    with open(conf_path, 'w') as f:
        f.write(conf_run_s)

#----------------
def mksh_job_run( conf_path, job_name, 
              nprc, elapse_time ):

  jobshell_header_s = batch_job_common.get_job_header(job_name, nprc, elapse_time)
  jobshell_s = f"""
  
SCALE_DG_INIT_BIN={SCALE_DG_BIN_PATH}/scale-dg_init
SCALE_DG_BIN={SCALE_DG_BIN_PATH}/scale-dg
llio_transfer ${{SCALE_DG_INIT_BIN}} ${{SCALE_DG_BIN}} *.conf

mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_INIT_BIN}} init.conf || exit 1

mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_BIN}} run.conf || exit 1  

llio_transfer --purge ${{SCALE_DG_INIT_BIN}} ${{SCALE_DG_BIN}} *.conf    
  """
  
  with open(conf_path, 'w') as f:
      f.write(jobshell_header_s + jobshell_s)

#----------------
def mksh_job_regrid( conf_path, job_name, 
              nprc, elapse_time ):

  jobshell_header_s = batch_job_common.get_job_header(job_name, nprc, elapse_time)
  jobshell_s = f"""
mkdir -p outdata/

SCALE_DG_REGRID_BIN={SCALE_DG_REGRID_BIN_PATH}/regrid_tool
llio_transfer ${{SCALE_DG_REGRID_BIN}} *.conf

mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_REGRID_BIN}} regrid.conf || exit 1

llio_transfer --purge 
llio_transfer ${{SCALE_DG_REGRID_BIN}} *.conf 
  """
  
  with open(conf_path, 'w') as f:
      f.write(jobshell_header_s + jobshell_s)
             
def mk_conf_jobsh( exp_name, exp_info, alph, exp_dir_postfix="" ):
  nprc = exp_info["nprc"]
  eh = exp_info["Eh"]
  ez = exp_info["Ez"]
  porder = exp_info["porder"]

  out_dir=f"./solid_body_rot_alph_{int(alph)}deg{exp_dir_postfix}/{ptracer_shape}/{exp_name}"

  print(out_dir)
  os.makedirs(out_dir, exist_ok=True)
  
  mkconf_init(f"{out_dir}/init.conf", 
              nprc, eh, ez, porder, 
              ptracer_shape, lonc, latc, rh, exp_info["initgp_porder"], 
              tau, alph )
  
  mkconf_run(f"{out_dir}/run.conf", 
              nprc, eh, ez, porder, 
              ptracer_shape, lonc, latc, rh, exp_info["initgp_porder"], 
              tau, alph,
              exp_info["dt"])        
  
  mkconf_regrid(f"{out_dir}/regrid.conf", 
              nprc, eh, ez, porder, 
              exp_info["regrid_nprcx"], exp_info["regrid_nprcy"], 
              exp_info["regrid_Ex"], exp_info["regrid_Ey"], exp_info["Ez"], exp_info["regrid_porder"] )
              
  mksh_job_run(f"{out_dir}/job_run.sh", f"ADV_H{eh}P{porder}", 
              nprc, exp_info["elapse_time"]) 
  
  mksh_job_regrid(f"{out_dir}/job_regrid.sh", f"ADVRG_H{eh}P{porder}", 
              exp_info["regrid_nprcx"]*exp_info["regrid_nprcy"], exp_info["regrid_elapse_time"])            
