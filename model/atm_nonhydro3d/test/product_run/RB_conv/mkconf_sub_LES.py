import os
import sys
import datetime
sys.path.append(os.path.join(os.path.dirname(__file__), './common'))
import batch_job_common

SCALE_DG_BIN_PATH="../../../"
SCALE_DG_REGRID_BIN_PATH="../../../../../../../../bin"

#----------------------

def mkconf_init( conf_path,
                nprcx, nprcy, eh, ez, porder, initgp_porder ): 
    conf_init_s = f"""#--- Configuration file for a test case of RB convection  -------
&PARAM_IO
 IO_LOG_BASENAME = 'init_LOG',
/
&PARAM_MKINIT
  initname = 'RB_conv', 
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
  CONST_OHM = 0.0D0, 
/
&PARAM_EXP
  THETA0 = 300.0D0, 
  DTHETA = -0.001D0, 
  x_c = 0.0D3, 
  y_c = 0.0D3,  
  z_c = 2.0D2,  
  r_x = 3.0D1, 
  r_y = 3.0D1,
  r_z = 3.0D1,
  BruntVaisalaFreq = 0.0D0, 
  IntrpPolyOrder_h = {initgp_porder}, 
  IntrpPolyOrder_v = {initgp_porder},   
/
&PARAM_CONST
  CONST_OHM = 0.0D0,
/
#** ATMOS ******************************************************
&PARAM_ATMOS
  ACTIVATE_FLAG = .true., 
  ATMOS_DYN_DO  = .true.
/
&PARAM_ATMOS_MESH
  NLocalMeshPerPrc = 1, 
  dom_xmin         = -1.6D3, 
  dom_xmax         =  1.6D3, 
  isPeriodicX      = .true.,
  dom_ymin         = -1.6D3,  
  dom_ymax         =  1.6D3,  
  isPeriodicY      = .true.,  
  dom_zmin         = 0.0D0, 
  dom_zmax         = 1.6D3,   
  NprcX            = {nprcx}, 
  NeX              = {eh},
  NprcY            = {nprcy},   
  NeY              = {eh},
  NeZ              = {ez},
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
                restart_in_basename, 
                start_month, start_day, start_hour, integ_sec, 
                dt_sec, dt_dyn_sec,                 
                nprcx, nprcy, eh, ez, porder, 
                modal_filter_flag, mf_alph, mf_ordh, mf_alpv, mf_ordv, 
                output_dtsec ): 
  
  if modal_filter_flag:
    mf_flag = ".true."
    param_mf = f"""
&PARAM_ATMOS_DYN_MODALFILTER
  MF_ETAC_h  = 0.0D0, 
  MF_ALPHA_h = {mf_alph}, 
  MF_ORDER_h = {mf_ordh},
  MF_ETAC_v  = 0.0D0, 
  MF_ALPHA_v = {mf_alpv}, 
  MF_ORDER_v = {mf_ordv},
/    
    """
  else:
    mf_flag = ".false."    
    param_mf = ""
    

  conf_run_s = f"""#--- Configuration file for a test case of RB convection  -------
&PARAM_RESTART
  IN_BASENAME  = "{restart_in_basename}",
  OUTPUT_FLAG  = .true., 
  OUT_BASENAME = 'restart'
/
&PARAM_TIME
  TIME_STARTDATE       = 1, {start_month}, {start_day}, {start_hour}, 0, 0,
  TIME_STARTMS         = 0.D0,
  TIME_DURATION        = {integ_sec},
  TIME_DURATION_UNIT   = 'SEC', 
  TIME_DT              = {dt_sec},
  TIME_DT_UNIT         = "SEC",
/
&PARAM_CONST
  CONST_OHM = 0.0D0, 
/
&PARAM_USER
  USER_do = .true., 
/
#** ATMOS ******************************************************
&PARAM_ATMOS
  ACTIVATE_FLAG       = .true., 
  TIME_DT             = {dt_sec}D0, 
  TIME_DT_UNIT        = 'SEC', 
  ATMOS_DYN_DO        = .true.
/
&PARAM_ATMOS_MESH
  NLocalMeshPerPrc = 1, 
  dom_xmin         = -1.6D3, 
  dom_xmax         =  1.6D3, 
  isPeriodicX      = .true.,
  dom_ymin         = -1.6D3,  
  dom_ymax         =  1.6D3,  
  isPeriodicY      = .true.,  
  dom_zmin         = 0.0D0, 
  dom_zmax         = 1.6D3,   
  NprcX            = {nprcx}, 
  NeX              = {eh},
  NprcY            = {nprcy},   
  NeY              = {eh},
  NeZ              = {ez},
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},  
/
&PARAM_ATMOS_VARS
  CHECK_RANGE = .true. ,
  CHECK_TOTAL = .false.,
/
#** ATMOS / DYN ******************************************************
&PARAM_ATMOS_DYN
  EQS_TYPE         = "NONHYDRO3D_RHOT_HEVE", 
  !-
  TINTEG_TYPE   = 'ERK_SSP_10s4o_2N', 
  TIME_DT          = {dt_dyn_sec}D0, 
  TIME_DT_UNIT     = 'SEC', 
  !-
  MODALFILTER_FLAG  = {mf_flag},
  NUMDIFF_FLAG      = .false., 
/
&PARAM_ATMOS_DYN_BND
  btm_vel_bc   = 'NOSLIP', 
  top_vel_bc   = 'NOSLIP', 
  btm_thermal_bc   = 'FIXVAL', 
  top_thermal_bc   = 'FIXVAL', 
  btm_thermal_fixval = 301.000D0,
  top_thermal_fixval = 284.382D0,    
  north_thermal_bc  = 'ADIABATIC',
  south_thermal_bc  = 'ADIABATIC', 
  west_thermal_bc   = 'ADIABATIC',
  east_thermal_bc   = 'ADIABATIC',
/
{param_mf}
#** ATMOS / PHYS / Turbulence ***********************************************
&PARAM_ATMOS_PHY_TB
  TIME_DT      = {dt_sec}D0,   
  TIME_DT_UNIT = 'SEC', 
  TB_TYPE      = 'SMAGORINSKY', 
/
&PARAM_ATMOS_PHY_TB_DGM_SMG
/
#*** OUTPUT *******************************************
&PARAM_FILE_HISTORY
 FILE_HISTORY_DEFAULT_BASENAME  = "history",
 FILE_HISTORY_DEFAULT_TINTERVAL = {output_dtsec}D0, 
 FILE_HISTORY_DEFAULT_TUNIT     = "SEC",
 FILE_HISTORY_DEFAULT_TAVERAGE  = .false.,
 FILE_HISTORY_DEFAULT_DATATYPE  = "REAL4",
 FILE_HISTORY_OUTPUT_STEP0      = .true.,
/
&HISTORY_ITEM name='DDENS'        /
&HISTORY_ITEM name='PT'           /
&HISTORY_ITEM name='U'            /
&HISTORY_ITEM name='V'            /
&HISTORY_ITEM name='W'            /
&HISTORY_ITEM name='T'            /

#*** Statistics *******************************************

&PARAM_MESHFIELD_STATISTICS
 use_globalcomm = .true.,
/
&PARAM_MONITOR 
  MONITOR_STEP_INTERVAL = 20,
/
&MONITOR_ITEM name='DDENS' /
&MONITOR_ITEM name='ENGT'  /
&MONITOR_ITEM name='ENGK'  /
&MONITOR_ITEM name='ENGI'  /
&MONITOR_ITEM name='ENGP'  /
    """
    
  with open(conf_path, 'w') as f:
      f.write(conf_run_s)

def mkconf_regrid( conf_path,
                nprcx, nex, nprcy, ney, nez, porder, 
                regrid_nprcx, regrid_nprcy, 
                regrid_nex, regrid_ney, regrid_nez, 
                regrid_porder ): 
    conf_run_s = f"""#--- Configuration file for a test case of RB convection  -------
&PARAM_IO
 IO_LOG_BASENAME = "regrid_LOG"
! IO_LOG_ALLNODE  = .true., 
/
&PARAM_CONST
  CONST_OHM = 0.0D0, 
/
&PARAM_REGRID_MESH
 in_MeshType  = "STRUCTURED3D", 
 out_MeshType = "STRUCTURED3D", 
/
&PARAM_REGRID_INTERP_FIELD
  !- input --------------------
  in_basename="history",      
  vars = "U", "V", "W", "DDENS", "PT", 
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
&PARAM_REGRID_INMESH3D_STRUCTURED
  NLocalMeshPerPrc = 1, 
  dom_xmin         = -1.6D3, 
  dom_xmax         =  1.6D3, 
  isPeriodicX      = .true.,
  dom_ymin         = -1.6D3,  
  dom_ymax         =  1.6D3,  
  isPeriodicY      = .true.,  
  dom_zmin         = 0.0D0, 
  dom_zmax         = 1.6D3,   
  NprcX            = {nprcx}, 
  NeX              = {nex},
  NprcY            = {nprcy},   
  NeY              = {ney},
  NeGZ             = {nez},
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},    
/
&PARAM_REGRID_OUTMESH3D_STRUCTURED
  NLocalMeshPerPrc = 1, 
  dom_xmin         = -1.6D3, 
  dom_xmax         =  1.6D3, 
  isPeriodicX      = .true.,
  dom_ymin         = -1.6D3,  
  dom_ymax         =  1.6D3,  
  isPeriodicY      = .true.,  
  dom_zmin         = 0.0D0, 
  dom_zmax         = 1.6D3, 
  NprcX       = {regrid_nprcx},       
  NeX         = {regrid_nex},           
  NprcY       = {regrid_nprcy}, 
  NeY         = {regrid_ney},    
  NeGZ        = {regrid_nez}, 
  PolyOrder_h = {regrid_porder}, 
  PolyOrder_v = {regrid_porder}, 
/
    """
    
    with open(conf_path, 'w') as f:
        f.write(conf_run_s)

#----------------

def mksh_job_run( conf_path, job_name, 
              nprc, elapse_time, init_flag ):

  if init_flag:
    init_cmd = f"""
mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_INIT_BIN}} init.conf || exit 1  
  """
  else:
    init_cmd = ""

  jobshell_header_s = batch_job_common.get_job_header(job_name, nprc, elapse_time)
  jobshell_s = f"""
SCALE_DG_INIT_BIN={SCALE_DG_BIN_PATH}/scale-dg_init
SCALE_DG_BIN={SCALE_DG_BIN_PATH}/scale-dg  
llio_transfer ${{SCALE_DG_INIT_BIN}} ${{SCALE_DG_BIN}} *.conf

{init_cmd}
mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_BIN}} run.conf || exit 1  

llio_transfer --purge ${{SCALE_DG_INIT_BIN}} ${{SCALE_DG_BIN}} *.conf    
  """
  
  with open(conf_path, 'w') as f:
      f.write(jobshell_header_s + jobshell_s)

#----------------
def mksh_job_regrid( conf_path, job_name, regrid_cnf, 
              nprc, elapse_time, outdir ):

  jobshell_header_s = batch_job_common.get_job_header(job_name, nprc, elapse_time)
  jobshell_s = f"""
mkdir -p {outdir}/
SCALE_DG_REGRID_BIN={SCALE_DG_REGRID_BIN_PATH}/regrid_tool
llio_transfer ${{SCALE_DG_REGRID_BIN}} *.conf

mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_REGRID_BIN}} {regrid_cnf} || exit 1

llio_transfer --purge ${{SCALE_DG_REGRID_BIN}} *.conf 
  """
  
  with open(conf_path, 'w') as f:
      f.write(jobshell_header_s + jobshell_s)
             
def mk_conf_sh( exp_name, exp_info, exp_postfix="" ):
    nprcx = exp_info["nprcx"]
    nprcy = exp_info["nprcy"]    
    eh = exp_info["Eh"]
    ez = exp_info["Ez"]
    porder = exp_info["porder"]

    out_dir_pref=f"./LES{exp_postfix}/{exp_name}"
    nprc = nprcx * nprcy
    
    print(out_dir_pref)    
    #-------
    os.makedirs(out_dir_pref+f"/run1", exist_ok=True)    
    mkconf_init(f"{out_dir_pref}/run1/init.conf", 
                nprcx, nprcy, eh, ez, porder, exp_info["initgp_porder"] )
    
    #---                     
    date_time0 = datetime.datetime(1, 1, 1, 0, 0, 0, 0) # + datetime.timedelta(days=exp_info["ini_day"])
    hr_per_run = exp_info["hr_per_run"]
    
    for runno in range(1,int( exp_info["integ_hour"] /hr_per_run)+1):
      date_time = date_time0 + datetime.timedelta(hours=hr_per_run*(runno-1))
      if runno > 1:
        restart_in_basename = f"../run{runno-1}/restart_{date_time.year:04}{date_time.month:02}{date_time.day:02}-{date_time.hour:02}0000.000"
      else:
        restart_in_basename = f"init_{date_time.year:04}{date_time.month:02}{date_time.day:02}-{date_time.hour:02}0000.000"
    
      os.makedirs(out_dir_pref+f"/run{runno}", exist_ok=True)
    
      mkconf_run(f"{out_dir_pref}/run{runno}/run.conf", 
            restart_in_basename, date_time.month, date_time.day, date_time.hour, 3600*hr_per_run, 
            exp_info["dt_sec"], exp_info["dt_dyn_sec"],            
            nprcx, nprcy, eh, ez, porder, 
            exp_info["modal_filter_flag"], exp_info["mf_alph"], exp_info["mf_ordh"], exp_info["mf_alpv"], exp_info["mf_ordv"], 
            exp_info["hist_int_sec"] )        
                    
      # mkconf_regrid(f"{out_dir_pref}/regrid.conf", 
      #                 nprc, eh, ez, porder,  
      #                 exp_info["regrid_nprcx"], exp_info["regrid_nprcy"], 
      #                 exp_info["regrid_Ex"], exp_info["regrid_Ey"], exp_info["Ez"], exp_info["regrid_porder"] )
                  
      mksh_job_run(f"{out_dir_pref}/run{runno}/job_run.sh", f"RB_L_P{porder}_rn{runno}", 
                  nprc, exp_info["elapse_time"], (runno==1) ) 
          
      # mksh_job_regrid(f"{out_dir_pref}/job_regrid.sh", f"REG_E{eh}P{porder}", "regrid.conf", 
      #                   exp_info["regrid_nprcx"]*exp_info["regrid_nprcy"], exp_info["regrid_elapse_time"], 
      #                   "outdata")
              
      # mksh_job_regrid(f"{out_dir_pref}/job_regrid_topo.sh", f"REGT_E{eh}P{porder}", "regrid_topo.conf", 
      #                   exp_info["regrid_nprcx"]*exp_info["regrid_nprcy"], exp_info["regrid_elapse_time"], 
      #                   "outdata")
  
#---------------------------------
