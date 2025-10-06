import os
import sys
import datetime
sys.path.append(os.path.join(os.path.dirname(__file__), './common'))
import batch_job_common

SCALE_DG_BIN_PATH="../../../"
SCALE_DG_REGRID_BIN_PATH="../../../../../../../../bin"

#----------------------

def mkconf_init( conf_path,
                nprcx, nprcy, ex, ey, ez, porder, initgp_porder, 
                ini_PT ): 
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
  TIME_STARTDATE       = 0001, 1, 1, 0, 0, 0,
  TIME_STARTMS         = 0.D0,
/
&PARAM_CONST
  CONST_OHM = 0.0D0, 
/
&PARAM_EXP
  U0     = 5D0, 
  THETA0 = {ini_PT}, 
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
  PERTURB_TYPE = "RANDOM", 
/
&PARAM_CONST
  CONST_OHM = 0.0D0,
/
#** ATMOS ******************************************************
&PARAM_ATMOS
  ACTIVATE_FLAG = .true., 
  ATMOS_DYN_DO  = .true., 
  ATMOS_PHY_TB_DO = .true.,     
/
&PARAM_ATMOS_MESH
  dom_xmin         = -1.6D3, 
  dom_xmax         =  1.6D3, 
  isPeriodicX      = .true.,
  dom_ymin         = -1.6D3,  
  dom_ymax         =  1.6D3,  
  isPeriodicY      = .true.,  
  dom_zmin         = 0.0D0, 
  dom_zmax         = 1.6D3,   
  NprcX            = {nprcx}, 
  NeX              = {ex},
  NprcY            = {nprcy},   
  NeY              = {ey},
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
                start_month, start_day, start_hour, start_min, integ_sec, 
                dt_sec, dt_dyn_sec,                 
                nprcx, nprcy, ex, ey, ez, porder, 
                modal_filter_flag, mf_alph, mf_ordh, mf_alpv, mf_ordv, 
                tb_filter_fac, 
                output_dtsec, 
                const_hflx,                
                stabCoef_bnd, 
                comp_comm_overlap ): 
 
  
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
    
  if comp_comm_overlap:
    eqs_type = "NONHYDRO3D_HEVE_CCO"
    hide_mpi_comm = "HIDE_MPI_COMM_FLAG = .true.,"
  else:
    eqs_type = "NONHYDRO3D_RHOT_HEVE"
    hide_mpi_comm = ""

  conf_run_s = f"""#--- Configuration file for a test case of RB convection  -------
&PARAM_RESTART
  IN_BASENAME  = "{restart_in_basename}",
  OUTPUT_FLAG  = .true., 
  OUT_BASENAME = 'restart'
/
&PARAM_TIME
  TIME_STARTDATE       = 1, {start_month}, {start_day}, {start_hour}, {start_min}, 0,
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
  BTM_BC_TYPE_HEAT = 'FixedFlux', 
  TOP_BC_TYPE_HEAT = 'FixedFlux', 
  BTM_FIXED_HEAT_FLUX = {const_hflx}, 
  TOP_FIXED_HEAT_FLUX = {const_hflx}, 
  StabCoef_bnd = {stabCoef_bnd}, 
  U0           = 5D0,   
/
#** ATMOS ******************************************************
&PARAM_ATMOS
  ACTIVATE_FLAG       = .true., 
  TIME_DT             = {dt_sec}, 
  TIME_DT_UNIT        = 'SEC', 
  ATMOS_DYN_DO        = .true., 
  ATMOS_PHY_TB_DO     = .true.,   
/
&PARAM_ATMOS_MESH
  dom_xmin         = -1.6D3, 
  dom_xmax         =  1.6D3, 
  isPeriodicX      = .true.,
  dom_ymin         = -1.6D3,  
  dom_ymax         =  1.6D3,  
  isPeriodicY      = .true.,  
  dom_zmin         = 0.0D0, 
  dom_zmax         = 1.6D3,   
  NprcX            = {nprcx}, 
  NeX              = {ex},
  NprcY            = {nprcy},   
  NeY              = {ey},
  NeZ              = {ez},
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},  
  COMM_USE_MPI_PC  = .true., 
  ELEMENT_OPERATION_TYPE = 'TensorProd3D',   
/
&PARAM_ATMOS_VARS
  CHECK_RANGE = .true. ,
  CHECK_TOTAL = .false.,
/
#** ATMOS / DYN ******************************************************
&PARAM_ATMOS_DYN
  EQS_TYPE         = '{eqs_type}', 
  !-
  TINTEG_TYPE   = 'ERK_SSP_10s4o_2N', 
  TIME_DT          = {dt_dyn_sec}, 
  TIME_DT_UNIT     = 'SEC', 
  !-
  MODALFILTER_FLAG  = {mf_flag},
  NUMDIFF_FLAG      = .false., 
  {hide_mpi_comm}
/
&PARAM_ATMOS_DYN_BND
  btm_vel_bc   = 'NOSLIP', 
  top_vel_bc   = 'NOSLIP', 
  btm_velx_fixval = 5.0D0, 
  top_velx_fixval = 5.0D0,   
  btm_thermal_bc   = 'ADIABATIC', 
  top_thermal_bc   = 'ADIABATIC', 
  north_thermal_bc  = 'ADIABATIC',
  south_thermal_bc  = 'ADIABATIC', 
  west_thermal_bc   = 'ADIABATIC',
  east_thermal_bc   = 'ADIABATIC',
/
{param_mf}
#** ATMOS / PHYS / Turbulence ***********************************************
&PARAM_ATMOS_PHY_TB
  TIME_DT      = {dt_sec},   
  TIME_DT_UNIT = 'SEC', 
  TB_TYPE      = 'SMAGORINSKY', 
/
&PARAM_ATMOS_PHY_TB_DGM_SMG
  filter_fac   = {tb_filter_fac}, 
/
#*** OUTPUT *******************************************
&PARAM_FILE_HISTORY
 FILE_HISTORY_DEFAULT_BASENAME  = "history",
 FILE_HISTORY_DEFAULT_TINTERVAL = {output_dtsec}, 
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

#----------------

def mksh_job_run( conf_path, job_name, 
              nprc, elapse_time, init_flag, node_shape ):

  if init_flag:
    init_cmd = f"""
mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_INIT_BIN}} init.conf || exit 1  
  """
  else:
    init_cmd = ""

  jobshell_header_s = batch_job_common.get_job_header(job_name, nprc, elapse_time, node_shape)
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
    ex = exp_info["Ex"]
    ey = exp_info["Ey"]    
    ez = exp_info["Ez"]
    porder = exp_info["porder"]
    runno_s = exp_info["runno_s"]       

    out_dir_pref=f"./LES{exp_postfix}/{exp_name}"
    nprc = nprcx * nprcy
        
    print(out_dir_pref)    
    #-------
    if runno_s == 1:
      os.makedirs(out_dir_pref+f"/run1", exist_ok=True)    
      mkconf_init(f"{out_dir_pref}/run1/init.conf", 
                  nprcx, nprcy, ex, ey, ez, porder, exp_info["initgp_porder"],
                  exp_info["ini_PT"] )
    
    #---                     
    date_time0 = datetime.datetime(1, 1, 1, 0, 0, 0, 0) + datetime.timedelta(seconds=exp_info["start_time_sec"])
    hr_per_run = exp_info["hr_per_run"]
    
    for runno in range(runno_s,runno_s+int( exp_info["integ_hour"] /hr_per_run)):
      date_time = date_time0 + datetime.timedelta(minutes=60*hr_per_run*(runno-runno_s))
      if runno > 1:
        restart_in_basename = f"../run{runno-1}/restart_{date_time.year:04}{date_time.month:02}{date_time.day:02}-{date_time.hour:02}{date_time.minute:02}00.000"
      else:
        restart_in_basename = f"init_{date_time.year:04}{date_time.month:02}{date_time.day:02}-{date_time.hour:02}{date_time.minute:02}00.000"
    
      os.makedirs(out_dir_pref+f"/run{runno}", exist_ok=True)
    
      mkconf_run(f"{out_dir_pref}/run{runno}/run.conf", 
            restart_in_basename, date_time.month, date_time.day, date_time.hour, date_time.minute, 3600*hr_per_run, 
            exp_info["dt_sec"], exp_info["dt_dyn_sec"],            
            nprcx, nprcy, ex, ey, ez, porder, 
            exp_info["modal_filter_flag"], exp_info["mf_alph"], exp_info["mf_ordh"], exp_info["mf_alpv"], exp_info["mf_ordv"], 
            exp_info["tb_filter_fac"], 
            exp_info["hist_int_sec"], 
            exp_info["const_hflx"],
            exp_info["StabCoef_bnd"], exp_info["comp_comm_overlap"] )
                    
                    
      if runno==1:
        batch_job_common.mkconf_regrid( 
                        f"{out_dir_pref}/run{runno_s}/regrid_bs.conf", 
                        nprcx, nprcy, ex, ey, ez, porder,  
                        exp_info["regrid_nprcx"], exp_info["regrid_nprcy"], 
                        exp_info["regrid_Ex"], exp_info["regrid_Ey"], ez, exp_info["regrid_porder"],
                        restart_in_basename, "'PRES_hyd', 'DENS_hyd'", "outdata/basic_state" )
                    
      batch_job_common.mkconf_regrid(f"{out_dir_pref}/run{runno}/regrid.conf", 
                      nprcx, nprcy, ex, ey, ez, porder,  
                      exp_info["regrid_nprcx"], exp_info["regrid_nprcy"], 
                      exp_info["regrid_Ex"], exp_info["regrid_Ey"], ez, exp_info["regrid_porder"],
                      "history", "'U', 'V', 'W', 'DDENS', 'PT'", "outdata/history" )
      
      mksh_job_run(f"{out_dir_pref}/run{runno}/job_run.sh", f"RB_L_P{porder}_rn{runno}", 
                  nprc, exp_info["elapse_time"], (runno==1), exp_info["node_shape"] ) 
          
      if runno==1:
        mksh_job_regrid(f"{out_dir_pref}/run{runno}/job_regrid_bs.sh", f"RG_L_P{porder}_rn{runno}", "regrid_bs.conf", 
                          exp_info["regrid_nprcx"]*exp_info["regrid_nprcy"], exp_info["regrid_elapse_time"], 
                          "outdata")
          
      mksh_job_regrid(f"{out_dir_pref}/run{runno}/job_regrid.sh", f"RG_L_P{porder}_rn{runno}", "regrid.conf", 
                        exp_info["regrid_nprcx"]*exp_info["regrid_nprcy"], exp_info["regrid_elapse_time"], 
                        "outdata")
  
#---------------------------------
