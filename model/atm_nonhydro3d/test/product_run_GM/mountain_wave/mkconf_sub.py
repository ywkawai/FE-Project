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
                Ueq, h0, ini_bg_force_flag, 
                highlat_tappering_flag ): 
  if ini_bg_force_flag:
    Ueq0 = "0D0"
  else:
    Ueq0 = Ueq
  
  if highlat_tappering_flag:
    param_highlat_tappering = f"""
  SCHAER_MERI_TAPER_FLAG = .true., 
  SCHAER_MERI_TAPER_TANH_Clat = 1.0471975511965976D0, ! 60 deg
  SCHAER_MERI_TAPER_TANH_LatWidth = 0.13962634015954636D0, ! 8 deg 
    """
  else:
    param_highlat_tappering = ""
    
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
  CONST_RADIUS = 38219.67606478705D0, ! R_earth / 166.7
/
&PARAM_EXP
  DCMIP_case = '2-1',  
  Ueq        = {Ueq}, 
  Ueq0       = {Ueq0}, 
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
  dom_zmax         = 30.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder}, 
  FZ               = {fz}, 
!  LumpedMassMatFlag = .true.,   
/
&PARAM_MKTOPO
  toponame     = 'SCHAER_GLOBAL',  
  OUT_BASENAME = 'TOPO', 
/
&PARAM_MKTOPO_SCHAER_GLOBAL
  SCHAER_Clon     = 3.141592653589793D0, ! PI
  SCHAER_Clat     = 0.D0, 
  SCHAER_R        = 5000.D0,
  SCHAER_LAMBDA   = 4000.D0, 
  SCHAER_SHAPE_ID = 1, 
  SCHAER_HEIGHT = {h0}.D0,  
  quasi_2D_flag = .true., 
{param_highlat_tappering}
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
                restart_in_basename, start_day, integ_hour, 
                nprc, neh, nez, porder, 
                fz, dt, mf_alph, mf_ordh, mf_alpv, mf_ordv,
                ini_bg_force_flag, ini_bg_force_tscale,
                ini_bg_force_turnoff_tstart, ini_bg_force_turnoff_tscale, 
                lsponge_tappering_flag ):                  
  
  if ini_bg_force_flag:
    param_ini_bg_force = f"""
  ini_bg_force_flag = .true., 
  ini_bg_force_tscale = {ini_bg_force_tscale}D0, 
  ini_bg_force_turnoff_tstart = {ini_bg_force_turnoff_tstart}D0, 
  ini_bg_force_turnoff_tscale = {ini_bg_force_turnoff_tscale}D0,   
    """
  else:
    param_ini_bg_force = ""
  
  if lsponge_tappering_flag:
    param_lsponge_tappering = f"""
  SL_MERI_TAPER_FLAG = .true., 
  SL_MERI_TAPER_TANH_Clat = 1.0471975511965976D0, ! 60 deg
  SL_MERI_TAPER_TANH_LatWidth = 0.13962634015954636D0, ! 8 deg    
    """
  else:
    param_lsponge_tappering = ""
    

  conf_run_s = f"""#--- Configuration file for a test case of mountain wave  -------
&PARAM_RESTART
  IN_BASENAME = "{restart_in_basename}",
  OUTPUT_FLAG = .true., 
  OUT_BASENAME = 'restart'
/
&PARAM_TIME
  TIME_STARTDATE       = 0000, 1, {start_day}, 0, 0, 0,
  TIME_STARTMS         = 0.D0,
  TIME_DURATION        = {integ_hour}, 
  TIME_DURATION_UNIT   = 'HOUR', 
  TIME_DT              = 0.375D0, 
  TIME_DT_UNIT         = 'SEC', 
/
&PARAM_CONST
  CONST_OHM = 0.0D0, 
  CONST_RADIUS = 38219.67606478705D0, ! R_earth / 166.7
/
&PARAM_EXP
  DCMIP_case = '2-1',  
  Ueq        = 20.0D0, 
/
&PARAM_USER
  USER_do = .true., 
  sponge_layer_flag = .true., 
  zTop              = 30D3, 
  SPONGE_HEIGHT     = 15D3, 
  SPONGE_EFOLD_SEC  = 100D0, 
  SPONGE_LAYER_FUNC_NAME    = "TANH",     
  lateral_sponge_layer_flag = .true., 
  LATERAL_SPONGE_EFOLD_SEC  = 200D0,   
{param_ini_bg_force}
{param_lsponge_tappering}
/
#** ATMOS ******************************************************
&PARAM_ATMOS
  ACTIVATE_FLAG       = .true., 
  ATMOS_MESH_TYPE     = 'GLOBAL',   
  TIME_DT             = 0.375D0, 
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
!  LumpedMassMatFlag = .true.,   
  TOPO_IN_BASENAME  = 'TOPO',   
/
&PARAM_ATMOS_VARS
  CHECK_RANGE = .true. ,
  CHECK_TOTAL = .false.,
/
#** ATMOS / DYN ******************************************************
&PARAM_ATMOS_DYN
  EQS_TYPE         = "GLOBALNONHYDRO3D_RHOT_HEVE", 
  !-
  TINTEG_TYPE      = 'ERK_SSP_10s4o_2N', ! [IMEX_ARK_232, IMEX_ARK324, ERK_SSP_3s3o]
  TIME_DT          = {dt}D0, 
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
&HISTORY_ITEM name='V'        /
&HISTORY_ITEM name='W'        /
&HISTORY_ITEM name='DDENS'    /
&HISTORY_ITEM name='THERM'    /
!&HISTORY_ITEM name='PRES'     /

#*** Statistics *******************************************

&PARAM_MESHFIELD_STATISTICS
 use_globalcomm = .true.,
/
&PARAM_MONITOR 
  MONITOR_STEP_INTERVAL = 30,
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
                regrid_porderh, regrid_porderv, regrid_fz, 
                output_dir ): 
    conf_run_s = f"""#--- Configuration file for a test case of mountain wave  -------
&PARAM_IO
 IO_LOG_BASENAME = "regrid_LOG"
! IO_LOG_ALLNODE  = .true., 
/
&PARAM_CONST
  CONST_OHM = 0.0D0, 
  CONST_RADIUS = 38219.67606478705D0, ! R_earth / 166.7
/
&PARAM_REGRID_MESH
 in_MeshType  = "CUBEDSPHERE3D", 
 out_MeshType = "LONLAT3D", 
/
&PARAM_REGRID_INTERP_FIELD
  !- input --------------------
  in_basename="history",      
  vars = "W", "U", "V", !"PRES_hyd", 
  !out_tinterval = 5,
/
&PARAM_REGRID_FILE
  !-- output ----------------
  out_basename="{output_dir}/history", 
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
  PolyOrder_h = {regrid_porderh}, 
  PolyOrder_v = {regrid_porderv}, 
  dom_xmin    =   0.0D0, 
  dom_xmax    = 360.0D0,   
  dom_ymin    = -90.0D0, 
  dom_ymax    =  90.0D0, 
  dom_zmin    = 0.0D0, 
  dom_zmax    = 30.0D3,   
  FZ          = {regrid_fz},  
/
!&PARAM_REGRID_VCOORD
!  vintrp_name     = 'HEIGHT', 
!  out_NeZ         = {regrid_nez},                 
!  out_PolyOrder_v = {porder},         
!  out_dom_vmin    = 0D0,         
!  out_dom_vmax    = 30D3, 
!  out_Fz = {fz}, 
!  in_topofile_basename = "outdata/topo", 
!  topo_varname         = "topo",           
!/
    """
    
    with open(conf_path, 'w') as f:
        f.write(conf_run_s)

def mkconf_regrid_topo( conf_path,
                nprc, neh, porder,  
                regrid_nprcx, regrid_nprcy, 
                regrid_nex, regrid_ney, 
                regrid_porder ): 
    conf_run_s = f"""#--- Configuration file for a test case of mountain wave  -------
&PARAM_IO
 IO_LOG_BASENAME = "regrid_topo_LOG"
! IO_LOG_ALLNODE  = .true., 
/
&PARAM_CONST
  CONST_OHM = 0.0D0, 
  CONST_RADIUS = 38219.67606478705D0, ! R_earth / 166.7
/
&PARAM_REGRID_MESH
 in_MeshType  = "CUBEDSPHERE2D", 
 out_MeshType = "LONLAT2D", 
/
&PARAM_REGRID_INTERP_FIELD
  !- input --------------------
  in_basename="./TOPO",     
  vars = "topo",  
/
&PARAM_REGRID_FILE
  !-- output ----------------
  out_basename="./outdata/topo", 
  out_UniformGrid=.false., 
/
&PARAM_REGRID_INMESH2D_CUBEDSPHERE
  NLocalMeshPerPrc = 1, 
  Nprc             = {nprc}, 
  NeGX             = {neh},
  NeGY             = {neh},
  PolyOrder_h      = {porder},
/
&PARAM_REGRID_OUTMESH2D_STRUCTURED
  NprcX       = {regrid_nprcx},       
  NeX         = {regrid_nex},           
  NprcY       = {regrid_nprcy}, 
  NeY         = {regrid_ney},    
  PolyOrder_h = {regrid_porder}, 
  dom_xmin    =   0.0D0, 
  dom_xmax    = 360.0D0,   
  dom_ymin    = -90.0D0, 
  dom_ymax    =  90.0D0, 
/
    """
    
    with open(conf_path, 'w') as f:
        f.write(conf_run_s)

#---------------- 
def mksh_job_run( conf_path, job_name, 
              nprc, elapse_time ):

  jobshell_header_s = batch_job_common.get_job_header(job_name, nprc, elapse_time)
  jobshell_s = f"""

SCALE_DG_INIT_BIN={SCALE_DG_BIN_PATH}/scale-dg_init_gm_mtwave_bc
SCALE_DG_BIN={SCALE_DG_BIN_PATH}/scale-dg_gm_mtwave_bc
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
  
mkdir -p {outdir}/

SCALE_DG_REGRID_BIN={SCALE_DG_REGRID_BIN_PATH}/regrid_tool
llio_transfer ${{SCALE_DG_REGRID_BIN}} *.conf

mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_REGRID_BIN}} {regrid_cnf} || exit 1

llio_transfer --purge ${{SCALE_DG_REGRID_BIN}} *.conf 
  """
  
  with open(conf_path, 'w') as f:
      f.write(jobshell_header_s + jobshell_s)
             
def mk_conf_sh( exp_name, exp_info ):
    nprc = exp_info["nprc"]
    eh = exp_info["Eh"]
    ez = exp_info["Ez"]
    fz = exp_info["fz"]
    porder = exp_info["porder"]
    h0 = exp_info["h0"]
    ini_bg_force_flag = exp_info["ini_bg_force_flag"]
    highlat_tappering_flag = exp_info["highlat_tappering_flag"]
    
    out_dir_pref=f"./rhot_heve/{exp_name}"

    print(out_dir_pref)
    os.makedirs(out_dir_pref, exist_ok=True)

    ueq = "20D0"
    mkconf_init(f"{out_dir_pref}/init.conf", 
                nprc, eh, ez, porder, 
                fz,
                ueq, h0, ini_bg_force_flag, 
                highlat_tappering_flag )
    
    mkconf_run(f"{out_dir_pref}/run.conf", 
               "init_00000101-000000.000", 1, exp_info["integ_hour"], 
                nprc, eh, ez, porder, 
                fz, exp_info["dt"], 
                exp_info["mf_alph"], exp_info["mf_ordh"], exp_info["mf_alpv"], exp_info["mf_ordv"], 
                ini_bg_force_flag, exp_info["ini_bg_force_tscale"], exp_info["ini_bg_force_turnoff_tstart"], exp_info["ini_bg_force_turnoff_tscale"], 
                highlat_tappering_flag )  
                    
    mkconf_regrid(f"{out_dir_pref}/regrid.conf", 
                    nprc, eh, ez, porder, fz, 
                    exp_info["regrid_nprcx"], exp_info["regrid_nprcy"], 
                    exp_info["regrid_Ex"], exp_info["regrid_Ey"], exp_info["Ez"], 
                    exp_info["regrid_porder"], exp_info["regrid_porder"], fz, 
                    "outdata" )
    
    mkconf_regrid(f"{out_dir_pref}/regrid_compari.conf", 
                    nprc, eh, ez, porder, fz, 
                    exp_info["regrid_nprcx_compari"], exp_info["regrid_nprcy_compari"], 
                    exp_info["regrid_Ex_compari"], exp_info["regrid_Ey_compari"], exp_info["regrid_Ez_compari"], 
                    exp_info["regrid_porderh_compari"], exp_info["regrid_porderv_compari"], exp_info["regrid_fz_compari"],
                    "outdata_compari" )

    mkconf_regrid_topo(f"{out_dir_pref}/regrid_topo.conf", 
                    nprc, eh, porder, 
                    exp_info["regrid_nprcx"], exp_info["regrid_nprcy"], 
                    exp_info["regrid_Ex"], exp_info["regrid_Ey"], exp_info["regrid_porder"] )
                
    mksh_job_run(f"{out_dir_pref}/job_run.sh", f"MTWV_Eh{eh}P{porder}", 
                      nprc, exp_info["elapse_time"]) 
        
    mksh_job_regrid(f"{out_dir_pref}/job_regrid.sh", f"REG_E{eh}P{porder}", "regrid.conf", 
                      exp_info["regrid_nprcx"]*exp_info["regrid_nprcy"], exp_info["regrid_elapse_time"], 
                      "outdata")
            
    mksh_job_regrid(f"{out_dir_pref}/job_regrid_topo.sh", f"REGT_E{eh}P{porder}", "regrid_topo.conf", 
                      exp_info["regrid_nprcx"]*exp_info["regrid_nprcy"], exp_info["regrid_elapse_time"], 
                      "outdata")

    mksh_job_regrid(f"{out_dir_pref}/job_regrid_compari.sh", f"REG_E{eh}P{porder}", "regrid_compari.conf", 
                      exp_info["regrid_nprcx_compari"]*exp_info["regrid_nprcy_compari"], exp_info["regrid_elapse_time"], 
                      "outdata_compari")
  
#---------------------------------
