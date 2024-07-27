import os
import math
import datetime
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../common'))
import batch_job_common
  
SCALE_DG_BIN_PATH="../../../"
SCALE_DG_REGRID_BIN_PATH="../../../../../../../../bin"

#----------------------

def mkconf_regrid_init( conf_path, 
                restart_in_basename, in_nprc, in_neh, in_nez, in_porder, in_fz, 
                restart_out_basename, out_nprc, out_neh, out_nez, out_porder, out_fz ): 
    conf_regrid_s = f"""#--- Configuration file for Held Suarez test  -------
&PARAM_IO
 IO_LOG_BASENAME = "regrid_LOG"
! IO_LOG_ALLNODE  = .true., 
/
&PARAM_REGRID_MESH
 in_MeshType  = "CUBEDSPHERE3D", 
 out_MeshType = "CUBEDSPHERE3D", 
/
&PARAM_REGRID_INTERP_FIELD
  !- input --------------------
  in_basename="{restart_in_basename}",      
  vars = "DDENS", "MOMX", "MOMY", "MOMZ", "THERM", "DENS_hyd", "PRES_hyd",  
/
&PARAM_REGRID_FILE
  !-- output ----------------
  out_basename="{restart_out_basename}", 
  out_UniformGrid=.false., 
/    
&PARAM_REGRID_INMESH3D_CUBEDSPHERE
  NLocalMeshPerPrc = 1, 
  Nprc             = {in_nprc}, 
  NeGX             = {in_neh},
  NeGY             = {in_neh},
  NeGZ             = {in_nez},
  dom_zmin         = 0.0D0, 
  dom_zmax         = 30.0D3, 
  PolyOrder_h      = {in_porder},
  PolyOrder_v      = {in_porder},
  Fz               = {in_fz}, 
/
&PARAM_REGRID_OUTMESH3D_CUBEDSPHERE
  NLocalMeshPerPrc = 1, 
  Nprc             = {out_nprc}, 
  NeGX             = {out_neh},
  NeGY             = {out_neh},
  NeGZ             = {out_nez},
  dom_zmin         = 0.0D0, 
  dom_zmax         = 30.0D3, 
  PolyOrder_h      = {out_porder},
  PolyOrder_v      = {out_porder},
  Fz               = {out_fz}, 
/   
    """
    
    with open(conf_path, 'w') as f:
        f.write(conf_regrid_s)

#----------------

def mkconf_run( conf_path, 
                restart_in_basename, start_year, start_mon, start_day, time_duration_day, 
                nprc, neh, nez, porder, 
                fz, dt, dt_dyn, mf_alph, mf_ordh, mf_alpv, mf_ordv,
                spinup1_flag ):
  
  if spinup1_flag:
    sponge_layer_param = f"""&PARAM_ATMOS_DYN_SPONGELAYER
  SL_WDAMP_TAU        = 28800.0D0, 
  SL_WDAMP_HEIGHT     = 12.0D3, 
  SL_HORIVELDAMP_FLAG = .true.
/
    """   
  else:
    sponge_layer_param = f"""&PARAM_ATMOS_DYN_SPONGELAYER
  SL_WDAMP_TAU        = 4320000.0D0, ! 50 Days
  SL_WDAMP_HEIGHT     = 12.0D3, 
  SL_HORIVELDAMP_FLAG = .false.
/   
    """   
    
  conf_run_s = f"""#--- Configuration file for Held-Suarez test  -------
&PARAM_RESTART
  IN_BASENAME = "{restart_in_basename}",
  OUTPUT_FLAG = .true., 
  OUT_BASENAME = 'restart'
/
&PARAM_TIME
  TIME_STARTDATE       = {start_year}, {start_mon}, {start_day}, 0, 0, 0,
  TIME_STARTMS         = 0.D0,
  TIME_DURATION        = {time_duration_day}D0, 
  TIME_DURATION_UNIT   = 'DAY', 
  TIME_DT              = {dt}D0, 
  TIME_DT_UNIT         = 'SEC', 
/
&PARAM_USER
  USER_do = .true., 
/
#** ATMOS ******************************************************
&PARAM_ATMOS
  ACTIVATE_FLAG       = .true., 
  ATMOS_MESH_TYPE     = 'GLOBAL',   
  TIME_DT             = {dt}D0, 
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
  dom_zmax         = 30.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},
  FZ               = {fz},   
  LumpedMassMatFlag = .true.,   
/
&PARAM_ATMOS_VARS
  CHECK_RANGE = .true. ,
  CHECK_TOTAL = .false.,
/
#** ATMOS / DYN ******************************************************
&PARAM_ATMOS_DYN
  EQS_TYPE         = "GLOBALNONHYDRO3D_RHOT_HEVI", 
  !-
  TINTEG_TYPE  = 'IMEX_ARK324', ! [IMEX_ARK_232, IMEX_ARK324, RK_TVD_3]
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
{sponge_layer_param}
#*** OUTPUT *******************************************
&PARAM_FILE_HISTORY
 FILE_HISTORY_DEFAULT_BASENAME  = "history",
 FILE_HISTORY_DEFAULT_TINTERVAL = 5D0, 
 FILE_HISTORY_DEFAULT_TUNIT     = "DAY",
 FILE_HISTORY_DEFAULT_TAVERAGE  = .false.,
 FILE_HISTORY_DEFAULT_DATATYPE  = "REAL4",
 FILE_HISTORY_OUTPUT_STEP0      = .true.,
/
&HISTORY_ITEM name='U'        /
&HISTORY_ITEM name='V'        /
&HISTORY_ITEM name='W'        /
&HISTORY_ITEM name='T'        /
&HISTORY_ITEM name='PRES'     /

#*** Statistics *******************************************

&PARAM_MESHFIELD_STATISTICS
 use_globalcomm = .true.,
/
&PARAM_MONITOR 
  MONITOR_STEP_INTERVAL = 6,
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
                regrid_porder ): 
  
  conf_run_s = f"""#--- Configuration file for a test case of Held Suarez test  -------
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
  vars = "W", "T", "PRES", 
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
  dom_zmax         = 30.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},
  FZ               = {fz},
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
  dom_zmax    = 30.0D3,   
  FZ               = {fz},  
/
    """
    
  with open(conf_path, 'w') as f:
    f.write(conf_run_s)

def mkconf_regrid_p( conf_path,
                nprc, neh, nez, porder, fz, 
                regrid_nprcx, regrid_nprcy, 
                regrid_nex, regrid_ney, regrid_nez, 
                regrid_porder, shallow_atm_approx, regrid_hori_uniform ): 
  
  if shallow_atm_approx:
    shallow_atm_approx_flag = ""
  else:
    shallow_atm_approx_flag = "SHALLOW_ATM_APPROX_FLAG = .false.,"
    
  if regrid_hori_uniform:
    out_UniformGrid_flag=".true."
    out_dir="./outdata_p_uniform"
  else:
    out_UniformGrid_flag=".false."    
    out_dir="./outdata_p"    
    
  conf_run_s = f"""#--- Configuration file for a test case of Held Suarez test  -------
&PARAM_IO
 IO_LOG_BASENAME = "regrid_p_LOG"
! IO_LOG_ALLNODE  = .true., 
/
&PARAM_REGRID_MESH
 in_MeshType  = "CUBEDSPHERE3D", 
 out_MeshType = "LONLAT3D", 
/
&PARAM_REGRID_INTERP_FIELD
  !- input --------------------
  in_basename="history",      
  vars = "W", "T", ! "PRES", 
  !out_tinterval = 5,
/
&PARAM_REGRID_FILE
  !-- output ----------------
  out_basename="{out_dir}/history", 
  out_UniformGrid={out_UniformGrid_flag}, 
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
  dom_zmax         = 30.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},
  FZ               = {fz},
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
  dom_zmax    = 30.0D3,   
  FZ               = {fz},  
/
&PARAM_REGRID_VCOORD
  vintrp_name     = 'PRESSURE', 
  out_NeZ         = 10,                
  out_PolyOrder_v = 3,         
  out_dom_vmin    = 1000D0,         
  out_dom_vmax    = 20D2,                  
  out_Fz          = 1000D2, 950D2, 850D2, 790D2, 680D2, 550D2, 400D2, 250D2, 100D2, 50D2, 30D2,    
  extrapolate     = .true.,
/
  """
    
  with open(conf_path, 'w') as f:
    f.write(conf_run_s)

def mkconf_regrid_spectra_p( conf_path,
                nprc, neh, nez, porder, fz ):

  out_dir="./outdata_p_spectra"    
     
  conf_run_s = f"""#--- Configuration file for a test case of Held Suarez test  -------
&PARAM_IO
 IO_LOG_BASENAME = "regrid_p_spectra_LOG"
! IO_LOG_ALLNODE  = .true., 
/
&PARAM_REGRID_MESH
 in_MeshType  = "CUBEDSPHERE3D", 
 out_MeshType = "CUBEDSPHERE3D",  
/
&PARAM_REGRID_INTERP_FIELD
  !- input --------------------
  in_basename="history",      
  vars = "U", "V", "W", "T", ! "PRES", 
  !out_tinterval = 5,
/
&PARAM_REGRID_FILE
  !-- output ----------------
  out_basename="{out_dir}/history", 
  out_UniformGrid=.false., 
/
&PARAM_REGRID_OPERATE_FIELD
  uvmet_conversion_flag = .false., 
/
&PARAM_REGRID_INMESH3D_CUBEDSPHERE
  NLocalMeshPerPrc = 1, 
  Nprc             = {nprc}, 
  NeGX             = {neh},
  NeGY             = {neh},
  NeGZ             = {nez},
  dom_zmin         = 0.0D0, 
  dom_zmax         = 30.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},
  FZ               = {fz},
/
&PARAM_REGRID_OUTMESH3D_CUBEDSPHERE
  NLocalMeshPerPrc = 1, 
  Nprc             = {nprc}, 
  NeGX             = {neh},
  NeGY             = {neh},
  NeGZ             = {nez},
  dom_zmin         = 0.0D0, 
  dom_zmax         = 30.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},
  FZ               = {fz},
/
&PARAM_REGRID_VCOORD
  vintrp_name     = 'PRESSURE', 
  out_NeZ         = 10,                
  out_PolyOrder_v = 3,         
  out_dom_vmin    = 1000D0,         
  out_dom_vmax    = 30D2,                  
  out_Fz          = 1000D2, 950D2, 850D2, 790D2, 680D2, 550D2, 400D2, 250D2, 100D2, 50D2, 30D2,    
  extrapolate     = .true.,
/
    """
    
  with open(conf_path, 'w') as f:
      f.write(conf_run_s)

def mkconf_sh_spectra( conf_path,
                nprc, neh, porder, Mt ):

  out_dir="./outdata_p_spectra"    
     
  conf_run_s = f"""#--- Configuration file for SH transform tools  -------
&PARAM_IO
 IO_LOG_BASENAME = "SH_LOG"
! IO_LOG_ALLNODE  = .true., 
/  
&PARAM_SH_TRANSFORM
  in_filebase="./outdata_p_spectra/history",  
  out_filebase="./outdata_p_spectra/spectral_data",  
  Mt                  = {Mt},
  LevelNum            = 3, 
  TARGET_LEVELS       = 850D2, 550D2, 250D2, 
  LEVEL_UNITS         = "hPa",
  VARS                = "U", "V", "W", "T", 
  TARGET_PROC_NUM_TOT = {nprc},
/
&PARAM_SH_MESH
  !---------------------
  Nprc             = {nprc}, 
  NeGX             = {neh},
  NeGY             = {neh},
  NeZ              = 10,
  dom_zmin         = 0.0D0, 
  dom_zmax         = 30.0D2,  
  PolyOrder_h      = {porder},
  PolyOrder_v      = 3, 
  IntrpPolyOrder_h = 11,   
  FZ               = 1000D2, 950D2, 850D2, 790D2, 680D2, 550D2, 400D2, 250D2, 100D2, 50D2, 30D2,  
/
    """
    
  with open(conf_path, 'w') as f:
      f.write(conf_run_s)

#----------------

def mksh_job_run( conf_path, job_name, 
              nprc, elapse_time,
                spinup1_flag ):

  jobshell_header_s = batch_job_common.get_job_header(job_name, nprc, elapse_time)
  
  if spinup1_flag:
    init_cmd = f"""mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_INIT_BIN}} init.conf || exit 1    
    """
  else:
    init_cmd = ""

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

  if len(outdir) > 0:
    mkdir_cmd = f"mkdir -p {outdir}"
  else:
    mkdir_cmd = ""

  jobshell_header_s = batch_job_common.get_job_header(job_name, nprc, elapse_time)
  jobshell_s = f"""
{mkdir_cmd}
SCALE_DG_REGRID_BIN={SCALE_DG_REGRID_BIN_PATH}/regrid_tool
llio_transfer ${{SCALE_DG_REGRID_BIN}} *.conf

mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_REGRID_BIN}} {regrid_cnf} || exit 1

llio_transfer --purge ${{SCALE_DG_REGRID_BIN}} *.conf 
  """
  
  with open(conf_path, 'w') as f:
      f.write(jobshell_header_s + jobshell_s)

#----------------
def mksh_job_spectra( conf_path, job_name, spectra_cnf, 
              nprc, elapse_time, outdir ):

  jobshell_header_s = batch_job_common.get_job_header(job_name, nprc, elapse_time)
  jobshell_s = f"""
mkdir -p {outdir}/
SCALE_DG_SPECTRA_BIN={SCALE_DG_REGRID_BIN_PATH}/sh_transform
llio_transfer ${{SCALE_DG_SPECTRA_BIN}} *.conf

mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_SPECTRA_BIN}} {spectra_cnf} || exit 1

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
    
    out_dir_pref=f"./rhot_hevi/{exp_name}"
    print(out_dir_pref)
    
    os.makedirs(out_dir_pref+"/init_regrid", exist_ok=True)
    mkconf_regrid_init(f"{out_dir_pref}/init_regrid/regrid_restart.conf", 
                exp_info["rg_in_basename"], exp_info["rg_in_nprc"], exp_info["rg_in_Eh"], exp_info["rg_in_Ez"], exp_info["rg_in_porder"], exp_info["rg_in_fz"],                 
                exp_info["rg_out_basename"], nprc, eh, ez, porder, fz )

    mksh_job_regrid(f"{out_dir_pref}/init_regrid/job_regrid_restart.sh", f"REGRS_E{eh}P{porder}", "regrid_restart.conf", 
                nprc, exp_info["regrid_elapse_time"], 
                "")
    
    #---                     
    date_time0 = datetime.datetime(1, 1, 1, 0, 0, 0, 0) + datetime.timedelta(days=exp_info["ini_day"])
    day_per_run = exp_info["day_per_run"]
    
    for runno in range(1,int( exp_info["integ_day"] /day_per_run)+1):
      date_time = date_time0 + datetime.timedelta(days=day_per_run*(runno-1))
      
      if runno > 1:
        prev_run_dir = f"../run{runno-1}"
      else:
        prev_run_dir = f"../init_regrid"
      
      os.makedirs(out_dir_pref+f"/run{runno}", exist_ok=True)
      mkconf_run(f"{out_dir_pref}/run{runno}/run.conf", 
                f"{prev_run_dir}/restart_{date_time.year:04}{date_time.month:02}{date_time.day:02}-000000.000", 
                date_time.year, date_time.month, date_time.day, day_per_run, 
                nprc, eh, ez, porder, fz, exp_info["dt"], exp_info["dt_dyn"],
                exp_info["mf_alph"], exp_info["mf_ordh"], exp_info["mf_alpv"], exp_info["mf_ordv"], False )
      
      mkconf_regrid(f"{out_dir_pref}/run{runno}/regrid.conf",
                  nprc, eh, ez, porder, fz, 
                  exp_info["regrid_nprcx"], exp_info["regrid_nprcy"], 
                  exp_info["regrid_Ex"], exp_info["regrid_Ey"], exp_info["Ez"], exp_info["regrid_porder"] )
      mkconf_regrid_p(f"{out_dir_pref}/run{runno}/regrid_p.conf",
                  nprc, eh, ez, porder, fz, 
                  exp_info["regrid_nprcx"], exp_info["regrid_nprcy"], 
                  exp_info["regrid_Ex"], exp_info["regrid_Ey"], exp_info["Ez"], exp_info["regrid_porder"], 
                  False, False )
      mkconf_regrid_p(f"{out_dir_pref}/run{runno}/regrid_p_uniform.conf",
                  nprc, eh, ez, porder, fz, 
                  exp_info["regrid_nprcx"], exp_info["regrid_nprcy"], 
                  exp_info["regrid_Ex"], exp_info["regrid_Ey"], exp_info["Ez"], exp_info["regrid_porder"], 
                  False, True )
      
      mkconf_regrid_spectra_p( f"{out_dir_pref}/run{runno}/regrid_p_spectra.conf",
                      nprc, eh, ez, porder, fz )
      mkconf_sh_spectra( f"{out_dir_pref}/run{runno}/spectral_analysis.conf",
                      nprc, eh, porder, exp_info["spectra_Mt"] )
            
      mksh_job_run(f"{out_dir_pref}/run{runno}/job_run.sh", f"HS_Eh{eh}P{porder}_sp1", 
                  nprc, exp_info["elapse_time"], False) 
    
      mksh_job_regrid(f"{out_dir_pref}/run{runno}/job_regrid.sh", f"REG_E{eh}P{porder}_{runno}", "regrid.conf", 
                  exp_info["regrid_nprcx"]*exp_info["regrid_nprcy"], exp_info["regrid_elapse_time"], 
                  "outdata")
      mksh_job_regrid(f"{out_dir_pref}/run{runno}/job_regrid_p.sh", f"REG_E{eh}P{porder}_p_{runno}", "regrid_p.conf", 
                  exp_info["regrid_nprcx"]*exp_info["regrid_nprcy"], exp_info["regrid_elapse_time"], 
                  "outdata_p")             
      mksh_job_regrid(f"{out_dir_pref}/run{runno}/job_regrid_p_uniform.sh", f"REGU_E{eh}P{porder}_p_{runno}", "regrid_p_uniform.conf", 
                  exp_info["regrid_nprcx"]*exp_info["regrid_nprcy"], exp_info["regrid_elapse_time"], 
                  "outdata_p_uniform")             
      mksh_job_regrid(f"{out_dir_pref}/run{runno}/job_regrid_p_spectra.sh", f"REGS_E{eh}P{porder}_p_{runno}", "regrid_p_spectra.conf", 
                  nprc, exp_info["regrid_elapse_time"], 
                  "outdata_p_spectra")             
      mksh_job_spectra(f"{out_dir_pref}/run{runno}/job_spectral_analysis.sh", f"SA_E{eh}P{porder}_p_{runno}", "spectral_analysis.conf", 
                  exp_info["spectra_nprc"], exp_info["spectra_elapse_time"], 
                  "outdata_p_spectra")             
  
#---------------------------------