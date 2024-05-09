import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../common'))
import batch_job_common

SCALE_DG_BIN_PATH="../../"
SCALE_DG_REGRID_BIN_PATH="../../../../../../../bin"

#----------------------

def mkconf_init( conf_path,
                nprc, neh, nez, porder, 
                fz,
                dom_xmin, dom_xmax, dom_dy,                 
                h0, topo_cx, 
                ini_bg_force_flag, 
                ini_gp_polyorder ): 
  
  if ini_bg_force_flag:
    ini_U0 = 0
  else:
    ini_U0 = 20
  
  conf_init_s = f"""#--- Configuration file for a test case of mountain wave  -------
&PARAM_IO
 IO_LOG_BASENAME = 'init_LOG',
/
&PARAM_MKINIT
  initname = 'mountain_wave', 
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
!  CONST_GRAV = 0.0D0,   
!  CONST_GRAV = 6.2D0,     
/
&PARAM_EXP
  U0               = {ini_U0}D0,  ! Wind speed at the equator
  TEMP0            = 300D0, 
  IntrpPolyOrder_h={ini_gp_polyorder}, IntrpPolyOrder_v={ini_gp_polyorder}, 
/
#** ATMOS ******************************************************
&PARAM_ATMOS
  ACTIVATE_FLAG = .true., 
  ATMOS_DYN_DO  = .true.
/
&PARAM_ATMOS_MESH
  dom_xmin         = {dom_xmin},  
  dom_xmax         = {dom_xmax}, 
  isPeriodicX      = .true.,
  dom_ymin         = 0.0D3,  
  dom_ymax         = {dom_dy},
  isPeriodicY      = .true.,   
  NprcX            = {nprc}, 
  NeX              = {neh},
  NprcY            = 1, 
  NeY              = 1, 
  NeZ              = {nez},   
  dom_zmin         = 0.0D0, 
  dom_zmax         = 30.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder}, 
  FZ               = {fz}, 
!  LumpedMassMatFlag = .true.,   
/
&PARAM_MKTOPO
  toponame     = 'SCHAER',  
  OUT_BASENAME = 'TOPO', 
/
&PARAM_MKTOPO_SCHAER
  SCHAER_CX     = {topo_cx}, 
  SCHAER_RX     = 5.0D3,
  SCHAER_LAMBDA = 4.0D3, 
  SCHAER_HEIGHT = {h0}D0,  
  IntrpPolyOrder_h ={ini_gp_polyorder},
/
#** ATMOS / DYN ******************************************************
&PARAM_ATMOS_DYN
  EQS_TYPE         = "NONHYDRO3D_HEVE", 
!  EQS_TYPE         = "NONHYDRO3D_ETOT_HEVE",   
  TINTEG_TYPE      = 'ERK_SSP_3s3o',    
/    
    """
    
  with open(conf_path, 'w') as f:
      f.write(conf_init_s)

#----------------

def mkconf_run( conf_path, 
                restart_in_basename, integ_hour, 
                nprc, neh, nez, porder, 
                fz,
                dom_xmin, dom_xmax, dom_dy,                 
                dt, dt_dyn, mf_alph, mf_ordh, mf_alpv, mf_ordv, 
                ini_bg_force_flag, ini_bg_force_tscale,
                ini_bg_force_turnoff_tstart, ini_bg_force_turnoff_tscale, 
                sl_apply_dens ): 
  
  if ini_bg_force_flag:
    param_ini_bg_force = f"""
  ini_bg_force_flag = .true., 
  ini_bg_force_tscale = {ini_bg_force_tscale}D0, 
  ini_bg_force_turnoff_tstart = {ini_bg_force_turnoff_tstart}D0, 
  ini_bg_force_turnoff_tscale = {ini_bg_force_turnoff_tscale}D0,   
    """
    if sl_apply_dens:
      param_ini_bg_force = f"""
  {param_ini_bg_force}
  SL_APPLY_DENS = .true., 
    """
      
  else:
    param_ini_bg_force = f"""
  ini_bg_force_flag = .false., 
    """
    
  conf_run_s = f"""#--- Configuration file for a test case of mountain wave  -------
&PARAM_RESTART
  IN_BASENAME = "{restart_in_basename}",
  OUTPUT_FLAG = .true., 
  OUT_BASENAME = 'restart'
/
&PARAM_TIME
  TIME_STARTDATE       = 0000, 1, 1, 0, 0, 0,
  TIME_STARTMS         = 0.D0,
  TIME_DURATION        = {integ_hour}, 
  TIME_DURATION_UNIT   = 'HOUR', 
  TIME_DT              = {dt}D0, 
  TIME_DT_UNIT         = 'SEC', 
/
&PARAM_CONST
  CONST_OHM = 0.0D0, 
!  CONST_GRAV = 0.0D0,     
!  CONST_GRAV = 6.2D0,     
/
&PARAM_USER
  USER_do = .true., 
/
&PARAM_USER_MTWAVE
  U0                = 20D0,   
  zTop              = 30D3, 
  SPONGE_HEIGHT     = 15D3, 
  SPONGE_EFOLD_SEC  = 100D0, 
  lateral_sponge_layer_flag = .true.,   
  SPONGE_LATERAL_WIDTH      = 120D3,   
  LATERAL_SPONGE_EFOLD_SEC  = 100D0,
{param_ini_bg_force}
/
#** ATMOS ******************************************************
&PARAM_ATMOS
  ACTIVATE_FLAG       = .true., 
  TIME_DT             = {dt}D0,  
  TIME_DT_UNIT        = 'SEC', 
  ATMOS_DYN_DO        = .true.
/
&PARAM_ATMOS_MESH
  dom_xmin         = {dom_xmin},  
  dom_xmax         = {dom_xmax}, 
  isPeriodicX      = .true.,
  dom_ymin         = 0.0D3,  
  dom_ymax         = {dom_dy},
  isPeriodicY      = .true.,   
  NprcX            = {nprc}, 
  NeX              = {neh},
  NprcY            = 1, 
  NeY              = 1, 
  NeZ              = {nez}, 
  dom_zmin         = 0.0D0, 
  dom_zmax         = 30.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder}, 
  FZ               = {fz},   
!  LumpedMassMatFlag = .true.,   
  TOPO_IN_BASENAME  = 'TOPO',   
/
&PARAM_ATMOS_VARS
  CHECK_RANGE = .true. ,
  CHECK_TOTAL = .false.,
/
#** ATMOS / DYN ******************************************************
&PARAM_ATMOS_DYN
  EQS_TYPE         = "NONHYDRO3D_RHOT_HEVE", 
!  EQS_TYPE         = "NONHYDRO3D_ETOT_HEVE",   
  !-
  TINTEG_TYPE      = 'ERK_SSP_10s4o_2N', ! [IMEX_ARK_232, IMEX_ARK324, ERK_SSP_3s3o]
  TIME_DT          = {dt_dyn}D0, 
  TIME_DT_UNIT     = 'SEC', 
  !-
  MODALFILTER_FLAG  = .true.,
  NUMDIFF_FLAG      = .false., 
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
#*** OUTPUT *******************************************
&PARAM_FILE_HISTORY
 FILE_HISTORY_DEFAULT_BASENAME  = "history",
 FILE_HISTORY_DEFAULT_TINTERVAL = 600D0, 
 FILE_HISTORY_DEFAULT_TUNIT     = "SEC",
 FILE_HISTORY_DEFAULT_TAVERAGE  = .false.,
 FILE_HISTORY_DEFAULT_DATATYPE  = "REAL8",
 FILE_HISTORY_OUTPUT_STEP0      = .true.,
/
&HISTORY_ITEM name='U'        /
&HISTORY_ITEM name='W'        /
&HISTORY_ITEM name='DDENS'    /
&HISTORY_ITEM name='THERM'    /
&HISTORY_ITEM name='GsqrtMOMX'    /
&HISTORY_ITEM name='GsqrtMOMZ'    /
&HISTORY_ITEM name='GsqrtMOMW'    /
&HISTORY_ITEM name='GsqrtDDENS'    /
!&HISTORY_ITEM name='GsqrtDPRES'    /
!&HISTORY_ITEM name='GsqrtG13DPRES'    /
!&HISTORY_ITEM name='DENS_dt_1'    /
!&HISTORY_ITEM name='DENS_dt_2'    /
!&HISTORY_ITEM name='DENS_dt_3'    /
!&HISTORY_ITEM name='DENS_dt_4'    /
!&HISTORY_ITEM name='DENS_dt_5'    /
!&HISTORY_ITEM name='DENS_dt_6'    /

#*** Statistics *******************************************

&PARAM_MESHFIELD_STATISTICS
 use_globalcomm = .true.,
/
&PARAM_MONITOR 
  MONITOR_STEP_INTERVAL = 100,
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
                regrid_nprc, regrid_neh, regrid_nez, regrid_porder, regrid_fz, 
                dom_xmin, dom_xmax, dom_dy, 
                outdata_dir
                ): 
    conf_run_s = f"""#--- Configuration file for a test case of mountain wave  -------
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
  vars = "U", "W", "THERM", "DDENS", "GsqrtDDENS", "GsqrtMOMX", "GsqrtMOMZ", "GsqrtMOMW", !"DENS_dt_1","DENS_dt_2","DENS_dt_3","DENS_dt_4","DENS_dt_5","DENS_dt_6",
  !out_tinterval = 5,
/
&PARAM_REGRID_FILE
  !-- output ----------------
  out_basename="{outdata_dir}/history", 
  out_dtype="REAL8",     
/
&PARAM_REGRID_INMESH3D_STRUCTURED
  dom_xmin         = {dom_xmin},  
  dom_xmax         = {dom_xmax}, 
  dom_ymin         = 0.0D3,  
  dom_ymax         = {dom_dy},
  NprcX            = {nprc}, 
  NeX              = {neh},
  NprcY            = 1, 
  NeY              = 1, 
  NeGZ             = {nez}, 
  dom_zmin         = 0.0D0, 
  dom_zmax         = 30.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder}, 
  FZ               = {fz},  
/
&PARAM_REGRID_OUTMESH3D_STRUCTURED
  dom_xmin         = {dom_xmin},  
  dom_xmax         = {dom_xmax}, 
  dom_ymin         = 0.0D3,  
  dom_ymax         = {dom_dy},
  NprcX            = {regrid_nprc}, 
  NeX              = {regrid_neh},
  NprcY            = 1, 
  NeY              = 1, 
  NeGZ             = {regrid_nez}, 
  dom_zmin         = 0.0D0, 
  dom_zmax         = 30.0D3, 
  PolyOrder_h      = {regrid_porder},
  PolyOrder_v      = {regrid_porder}, 
  FZ               = {regrid_fz},          
/
    """
    
    with open(conf_path, 'w') as f:
        f.write(conf_run_s)

#---------------- 
def mksh_job_run( conf_path, job_name, 
              nprc, elapse_time ):

  jobshell_header_s = batch_job_common.get_job_header(job_name, nprc, elapse_time)
  jobshell_s = f"""

SCALE_DG_INIT_BIN={SCALE_DG_BIN_PATH}/scale-dg_init_rm_mtwave
SCALE_DG_BIN={SCALE_DG_BIN_PATH}/scale-dg_rm_mtwave
SCALE_DG_REGRID_BIN={SCALE_DG_REGRID_BIN_PATH}/regrid_tool

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
  
rm -rf {outdir}/
mkdir -p {outdir}/

SCALE_DG_REGRID_BIN={SCALE_DG_REGRID_BIN_PATH}/regrid_tool

llio_transfer ${{SCALE_DG_REGRID_BIN}} *.conf

mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_REGRID_BIN}} {regrid_cnf} || exit 1

llio_transfer --purge ${{SCALE_DG_REGRID_BIN}} *.conf 
  """
  
  with open(conf_path, 'w') as f:
      f.write(jobshell_header_s + jobshell_s)
             
def mk_conf_sh( exp_name, exp_info, out_dir_pref0="./rhot_heve_v2" ):
    nprc = exp_info["nprc"]
    eh = exp_info["Eh"]
    ez = exp_info["Ez"]
    fz = exp_info["fz"]
    dom_xmin = exp_info["dom_xmin"]
    dom_xmax = exp_info["dom_xmax"]
    dom_dy = exp_info["dom_dy"]    
    porder = exp_info["porder"]

    out_dir_pref=f"{out_dir_pref0}/{exp_name}"

    print(out_dir_pref)
    os.makedirs(out_dir_pref, exist_ok=True)
    
    mkconf_init(f"{out_dir_pref}/init.conf", 
                nprc, eh, ez, porder, fz, 
                dom_xmin, dom_xmax, dom_dy, 
                exp_info["h0"], exp_info["topo_cx"], 
                exp_info["ini_bg_force_flag"],
                exp_info["ini_gp_polyorder"] )
    
    mkconf_run(f"{out_dir_pref}/run.conf", 
               "init_00000101-000000.000", exp_info["integ_hour"], 
                nprc, eh, ez, porder, fz, 
                dom_xmin, dom_xmax, dom_dy,                 
                exp_info["dt"], exp_info["dt_dyn"], exp_info["mf_alph"], exp_info["mf_ordh"], exp_info["mf_alpv"], exp_info["mf_ordv"], 
                exp_info["ini_bg_force_flag"], exp_info["ini_bg_force_tscale"], exp_info["ini_bg_force_turnoff_tstart"], exp_info["ini_bg_force_turnoff_tscale"],
                exp_info["sponge_layer_apply_dens_flag"] )  
                    
    mkconf_regrid(f"{out_dir_pref}/regrid_compari.conf", 
                    nprc, eh, ez, porder, fz, 
                    exp_info["regrid_nprc"], exp_info["regrid_Ex"], exp_info["regrid_Ez"], exp_info["regrid_porder"], exp_info["regrid_fz"], 
                    dom_xmin, dom_xmax, exp_info["regrid_dom_dy"],                 
                    "outdata_compari" )
                
    mksh_job_run(f"{out_dir_pref}/job_run.sh", f"MTWV_Eh{nprc*eh}P{porder}", 
                      nprc, exp_info["elapse_time"]) 
        
    mksh_job_regrid(f"{out_dir_pref}/job_regrid_compari.sh", f"REG_E{nprc*eh}P{porder}", "regrid_compari.conf", 
                      exp_info["regrid_nprc"], exp_info["regrid_elapse_time"], 
                      "outdata_compari")
  
#---------------------------------
