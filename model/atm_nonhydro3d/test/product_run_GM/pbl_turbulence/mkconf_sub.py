import os
import math
import datetime
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../common'))
import batch_job_common

SCALE_DG_BIN_PATH="../../../"
SCALE_DG_REGRID_BIN_PATH="../../../../../../../../bin"

#----------------------

def mkconf_init( conf_path,
                nprc, neh, nez, porder, fz, rplanet, shallow_atm_approx ): 
  if shallow_atm_approx:
    shallow_atm_approx_flag = ""
  else:
    shallow_atm_approx_flag = "SHALLOW_ATM_APPROX_FLAG = .false.,"
  
  conf_init_s = f"""#--- Configuration file for a test case of PBL turbulence  -------
&PARAM_IO
 IO_LOG_BASENAME = 'init_LOG',
/
&PARAM_MKINIT
  initname = 'pbl_turblence', 
/
&PARAM_RESTART
  OUTPUT_FLAG = .true., 
  OUT_BASENAME = 'init'
/
&PARAM_TIME
  TIME_STARTDATE       = 0000, 1, 1, 0, 0, 0,
  TIME_STARTMS         = 0.D0,
/
&PARAM_EXP
  ENV_U          = 0D0,
  ENV_THETA_SFC  = 299.5D0, 
  ENV_THETA_LAPS = 4.0D-3, 
  RANDOM_THETA   = 1.0D0, 
/
&PARAM_CONST
  CONST_OHM = 0.0D0, 
  CONST_RADIUS = {rplanet}D0,  
/
&PARAM_CONST
  CONST_OHM = 0.0D0,
/
#** ATMOS ******************************************************
&PARAM_ATMOS
  ATMOS_MESH_TYPE = 'GLOBAL', 
  ACTIVATE_FLAG = .true., 
  ATMOS_DYN_DO  = .true.
/
&PARAM_ATMOS_MESH
  {shallow_atm_approx_flag}
  NLocalMeshPerPrc = 1, 
  Nprc             = {nprc}, 
  NeGX             = {neh},
  NeGY             = {neh},
  NeZ              = {nez},
  dom_zmin         = 0.0D0, 
  dom_zmax         = 3.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder}, 
  Fz               = {fz}, 
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
                restart_in_basename, start_hour, start_min, integ_hour, 
                nprc, neh, nez, porder, 
                dt, dt_dyn, 
                mf_alph, mf_ordh, mf_alpv, mf_ordv, 
                fz, rplanet, hist_int_sec, shallow_atm_approx ): 
  
  if shallow_atm_approx:
    shallow_atm_approx_flag = ""
  else:
    shallow_atm_approx_flag = "SHALLOW_ATM_APPROX_FLAG = .false.,"

  conf_run_s = f"""#--- Configuration file for a test case of PBL turbulence  -------
&PARAM_RESTART
  IN_BASENAME  = "{restart_in_basename}",
  OUTPUT_FLAG  = .true., 
  OUT_BASENAME = 'restart'
/
&PARAM_TIME
  TIME_STARTDATE       = 0000, 1, 1, {start_hour}, {start_min}, 0,
  TIME_STARTMS         = 0.D0,
  TIME_DURATION        = {integ_hour}D0, 
  TIME_DURATION_UNIT   = 'HOUR', 
  TIME_DT              = {dt}D0, 
  TIME_DT_UNIT         = 'SEC', 
/
&PARAM_CONST
  CONST_OHM = 0.0D0, 
  CONST_RADIUS = {rplanet}, 
/
&PARAM_USER
  USER_do = .true., 
/
#** ATMOS ******************************************************
&PARAM_ATMOS
  ACTIVATE_FLAG        = .true., 
  ATMOS_MESH_TYPE      = 'GLOBAL',   
  TIME_DT              = {dt}D0, 
  TIME_DT_UNIT         = 'SEC', 
  TIME_DT_RESTART      = 900D0, 
  TIME_DT_RESTART_UNIT = "SEC",      
  ATMOS_DYN_DO         = .true., 
  ATMOS_PHY_SF_DO      = .true.,
  ATMOS_PHY_TB_DO      = .true.,   
/
&PARAM_ATMOS_MESH
  {shallow_atm_approx_flag}
  NLocalMeshPerPrc = 1, 
  Nprc             = {nprc}, 
  NeGX             = {neh},
  NeGY             = {neh},
  NeZ              = {nez},
  dom_zmin         = 0.0D0, 
  dom_zmax         = 3.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},
  FZ               = {fz}, 
  LumpedMassMatFlag = .false.,   
/
&PARAM_ATMOS_VARS
  CHECK_RANGE = .true. ,
  CHECK_TOTAL = .false.,
/
#** ATMOS / DYN ******************************************************
&PARAM_ATMOS_DYN
  EQS_TYPE         = "GLOBALNONHYDRO3D_RHOT_HEVE", 
  !-
  TINTEG_TYPE      = 'ERK_SSP_10s4o_2N',
  TIME_DT          = {dt_dyn}D0, 
  TIME_DT_UNIT     = 'SEC', 
  !-
  MODALFILTER_FLAG  = .true.,
  NUMDIFF_FLAG      = .false., 
  SPONGELAYER_FLAG  = .true.,   
/
&PARAM_ATMOS_DYN_BND
  btm_vel_bc   = 'SLIP', 
  top_vel_bc   = 'SLIP', 
  btm_thermal_bc   = 'ADIABATIC', 
  top_thermal_bc   = 'ADIABATIC', 
/
&PARAM_ATMOS_DYN_MODALFILTER
  MF_ETAC_h  = 0.0D0, 
  MF_ALPHA_h = {mf_alph}, 
  MF_ORDER_h = {mf_ordh},
  MF_ETAC_v  = 0.0D0, 
  MF_ALPHA_v = {mf_alpv}, 
  MF_ORDER_v = {mf_ordv},
/
&PARAM_ATMOS_DYN_SPONGELAYER
  SL_WDAMP_TAU     = 10D0,
  SL_WDAMP_HEIGHT  = 2000D0, 
/
#** ATMOS / PHYS / SFC ******************************************************
&PARAM_ATMOS_PHY_SFC
  TIME_DT       = {dt}D0, 
  TIME_DT_UNIT  = 'SEC', 
  SFCFLX_TYPE   = 'CONST'
/
&PARAM_ATMOS_PHY_SF_CONST
 ATMOS_PHY_SF_Const_SH =   200.D0,
 ATMOS_PHY_SF_Const_LH =     0.D0,
 ATMOS_PHY_SF_Const_Cm = 0.0011D0,
/
#** ATMOS / PHYS / Turbulence ***********************************************
&PARAM_ATMOS_PHY_TB
  TIME_DT      = {dt}D0, 
  TIME_DT_UNIT = 'SEC', 
  TB_TYPE      = 'SMAGORINSKY_GLOBAL'
/
#*** OUTPUT *******************************************
&PARAM_FILE_HISTORY
 FILE_HISTORY_DEFAULT_BASENAME  = "history",
 FILE_HISTORY_DEFAULT_TINTERVAL = {hist_int_sec}D0, 
 FILE_HISTORY_DEFAULT_TUNIT     = "SEC",
 FILE_HISTORY_DEFAULT_TAVERAGE  = .false.,
 FILE_HISTORY_DEFAULT_DATATYPE  = "REAL4",
 FILE_HISTORY_OUTPUT_STEP0      = .true.,
/
&HISTORY_ITEM name='DDENS'        /
&HISTORY_ITEM name='U'            /
&HISTORY_ITEM name='V'            /
&HISTORY_ITEM name='W'            /
&HISTORY_ITEM name='PT'           /

#*** Statistics *******************************************

&PARAM_MESHFIELD_STATISTICS
 use_globalcomm = .true.,
/
&PARAM_MONITOR 
  MONITOR_STEP_INTERVAL = 50,
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
                nprc, neh, nez, porder, fz, 
                regrid_nprcx, regrid_nprcy, 
                regrid_nex, regrid_ney, regrid_nez, 
                regrid_porder, rplanet, 
                out_dir, uniform_grid_flag, shallow_atm_approx ): 
  
  if shallow_atm_approx:
    shallow_atm_approx_flag = ""
  else:
    shallow_atm_approx_flag = "SHALLOW_ATM_APPROX_FLAG = .false.,"
  
  conf_run_s = f"""#--- Configuration file for a test case of PBL turbulence  -------
&PARAM_IO
 IO_LOG_BASENAME = "regrid_LOG"
! IO_LOG_ALLNODE  = .true., 
/
&PARAM_CONST
  CONST_OHM = 0.0D0, 
  CONST_RADIUS = {rplanet}, ! R_earth / 166.7
/
&PARAM_REGRID_MESH
 in_MeshType  = "CUBEDSPHERE3D", 
 out_MeshType = "LONLAT3D", 
/
&PARAM_REGRID_INTERP_FIELD
  !- input --------------------
  in_basename="history",      
  vars = "W", "U", "V", "DDENS", "PT", 
  !out_tinterval = 5,
/
&PARAM_REGRID_FILE
  !-- output ----------------
  out_basename="./{out_dir}/history", 
  out_UniformGrid={uniform_grid_flag}, 
/
&PARAM_REGRID_OPERATE_FIELD
  uvmet_conversion_flag = .true., 
/
&PARAM_REGRID_INMESH3D_CUBEDSPHERE
  {shallow_atm_approx_flag}
  NLocalMeshPerPrc = 1, 
  Nprc             = {nprc}, 
  NeGX             = {neh},
  NeGY             = {neh},
  NeGZ             = {nez},
  dom_zmin         = 0.0D0, 
  dom_zmax         = 3.0D3, 
  FZ               = {fz}, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},
/
&PARAM_REGRID_OUTMESH3D_STRUCTURED
  {shallow_atm_approx_flag}
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
  dom_zmax    = 3.0D3,   
  FZ          = {fz},   
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

def mksh_job_spectra_analysis( conf_path, job_name, analysis_cnf, 
              nprc, elapse_time, outdir ):

  jobshell_header_s = batch_job_common.get_job_header(job_name, nprc, elapse_time)
  jobshell_s = f"""
mkdir -p {outdir}/
SCALE_DG_SPECTRA_BIN={SCALE_DG_REGRID_BIN_PATH}/sh_transform
llio_transfer ${{SCALE_DG_SPECTRA_BIN}} *.conf

mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_SPECTRA_BIN}} {analysis_cnf} || exit 1

llio_transfer --purge ${{SCALE_DG_SPECTRA_BIN}} *.conf 
  """
  
  with open(conf_path, 'w') as f:
      f.write(jobshell_header_s + jobshell_s)
             
def mk_conf_sh( exp_name, exp_info ):
    nprc = exp_info["nprc"]
    eh = exp_info["Eh"]
    ez = exp_info["Ez"]
    fz = exp_info["fz"]    
    porder = exp_info["porder"]
    rplanet = exp_info["rplanet"]
    run_num = exp_info["run_num"]
    shallow_atm_approx = exp_info["shallow_atm_approx"] 

    rplanet_km = '{:.1f}'.format(rplanet/1000)
    out_dir_pref0=f"./rp{rplanet_km}km/{exp_name}"

    print(out_dir_pref0)
    os.makedirs(out_dir_pref0, exist_ok=True)
        
    date_time0 = datetime.datetime(1, 1, 1, 0, 0, 0, 0)
    hr_per_run = exp_info["hour_per_run"]
    
    for runno in range(1,run_num+1):
      date_time = date_time0 + datetime.timedelta(hours=hr_per_run*(runno-1))
      if runno >= exp_info["runno_analysis"]:
        hist_int_sec = exp_info["hist_int_sec_analysis"] 
      else:
        hist_int_sec = exp_info["hist_int_sec"]
      
      out_dir_pref = f"{out_dir_pref0}/run{runno}/"
      os.makedirs(out_dir_pref, exist_ok=True)
      
      if runno ==1:
        restart_in_basename = "init_00000101-000000.000"
        mkconf_init(f"{out_dir_pref0}/run1/init.conf", 
                    nprc, eh, ez, porder, fz, rplanet, shallow_atm_approx ) 
        mf_alph =  exp_info["mf_alph_ini"]; mf_alpv =  exp_info["mf_alpv_ini"]
      else:
        restart_in_basename = f"../run{runno-1}/restart_0000{date_time.month:02}{date_time.day:02}-{date_time.hour:02}{date_time.minute:02}00.000" 
        mf_alph =  exp_info["mf_alph"]; mf_alpv =  exp_info["mf_alpv"]

      mkconf_run(f"{out_dir_pref}/run.conf", 
                restart_in_basename, date_time.hour, date_time.minute, hr_per_run, 
                nprc, eh, ez, porder, 
                exp_info["dt"], exp_info["dt_dyn"],
                mf_alph, exp_info["mf_ordh"], mf_alpv, exp_info["mf_ordv"], 
                fz, rplanet, hist_int_sec, shallow_atm_approx )        
                      
      mkconf_regrid(f"{out_dir_pref}/regrid.conf", 
                      nprc, eh, ez, porder, fz, 
                      exp_info["regrid_nprcx"], exp_info["regrid_nprcy"], 
                      exp_info["regrid_Ex"], exp_info["regrid_Ey"], exp_info["Ez"], exp_info["regrid_porder"],
                      rplanet, "./outdata/", ".false.", shallow_atm_approx )

      # mkconf_regrid(f"{out_dir_pref}/regrid_uniform.conf", 
      #                 nprc, eh, ez, porder, fz, 
      #                 exp_info["regrid_nprcx"], exp_info["regrid_nprcy"], 
      #                 exp_info["regrid_Ex"], exp_info["regrid_Ey"], exp_info["Ez"], exp_info["regrid_porder"],
      #                 rplanet, "./output_uniform/", ".true.", shallow_atm_approx )
      
                  
      mksh_job_run(f"{out_dir_pref}/job_run.sh", f"PBL_Eh{eh}P{porder}", 
                        nprc, exp_info["elapse_time"]) 
          
      mksh_job_regrid(f"{out_dir_pref}/job_regrid.sh", f"REG_E{eh}P{porder}", "regrid.conf", 
                        exp_info["regrid_nprcx"]*exp_info["regrid_nprcy"], exp_info["regrid_elapse_time"], 
                        "outdata")
              
      # mksh_job_regrid(f"{out_dir_pref}/job_regrid_uniform.sh", f"REG_E{eh}P{porder}", "regrid_uniform.conf", 
      #                   exp_info["regrid_nprcx"]*exp_info["regrid_nprcy"], exp_info["regrid_elapse_time"], 
      #                   "outdata")
  
#---------------------------------
