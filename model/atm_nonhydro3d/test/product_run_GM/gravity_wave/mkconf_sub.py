import os
import math

SCALE_DG_BIN_PATH="../../"
SCALE_DG_REGRID_BIN_PATH="../../../../../../../bin"

#----------------------

def mkconf_init( conf_path,
                nprc, neh, nez, porder, initgp_porder ): 
    conf_init_s = f"""#--- Configuration file for a test case of gravity wave  -------
&PARAM_IO
 IO_LOG_BASENAME = 'init_LOG',
/
&PARAM_MKINIT
  initname = 'inertia_gravity_wave', 
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
  BruntVaisalaFreq = 1.0D-2, 
  THETA0           = 300D0, 
  InitAtmType      = 'ISOTHERMAL', 
  InitPerturbType  = 'GAUSSIAN',
  rh               = 1274244D0, ! RPlanet / 5.0
  DTHETA           = 1D-3, 
  nv               = 1,
  Ini_GP_PolyOrder_h = {initgp_porder}, 
  Ini_GP_PolyOrder_v = {initgp_porder}, 
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
                restart_in_basename, start_day, 
                nprc, neh, nez, porder, 
                dt, dt_dyn, 
                modal_filter_flag, mf_alph, mf_ordh, mf_alpv, mf_ordv ): 
  
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
    

  conf_run_s = f"""#--- Configuration file for a test case of gravity wave  -------
&PARAM_RESTART
  IN_BASENAME  = "{restart_in_basename}",
  OUTPUT_FLAG  = .true., 
  OUT_BASENAME = 'restart'
/
&PARAM_TIME
  TIME_STARTDATE       = 0000, 1, {start_day}, 0, 0, 0,
  TIME_STARTMS         = 0.D0,
  TIME_DURATION        = 2.0D0, 
  TIME_DURATION_UNIT   = 'DAY', 
  TIME_DT              = {dt}D0, 
  TIME_DT_UNIT         = 'SEC', 
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
  dom_zmax         = 10.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},
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
  TINTEG_TYPE      = 'IMEX_ARK324', ! [IMEX_ARK_232, IMEX_ARK324, ERK_SSP_3s3o]
  TIME_DT          = {dt_dyn}D0, 
  TIME_DT_UNIT     = 'SEC', 
  !-
  MODALFILTER_FLAG  = {mf_flag},
  NUMDIFF_FLAG      = .false., 
/
&PARAM_ATMOS_DYN_BND
  btm_vel_bc   = 'SLIP', 
  top_vel_bc   = 'SLIP', 
  btm_thermal_bc   = 'ADIABATIC', 
  top_thermal_bc   = 'ADIABATIC', 
/
{param_mf}
#*** OUTPUT *******************************************
&PARAM_FILE_HISTORY
 FILE_HISTORY_DEFAULT_BASENAME  = "history",
 FILE_HISTORY_DEFAULT_TINTERVAL = 14400D0, 
 FILE_HISTORY_DEFAULT_TUNIT     = "SEC",
 FILE_HISTORY_DEFAULT_TAVERAGE  = .false.,
 FILE_HISTORY_DEFAULT_DATATYPE  = "REAL8",
 FILE_HISTORY_OUTPUT_STEP0      = .true.,
/
&HISTORY_ITEM name='U'        /
&HISTORY_ITEM name='V'        /
&HISTORY_ITEM name='W'        /
!&HISTORY_ITEM name='PT_diff'  /
&HISTORY_ITEM name='DDENS'    /
&HISTORY_ITEM name='THERM'    /

#*** Statistics *******************************************

&PARAM_MESHFIELD_STATISTICS
 use_globalcomm = .true.,
/
&PARAM_MONITOR 
  MONITOR_STEP_INTERVAL = 2,
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
                nprc, neh, nez, porder, 
                regrid_nprcx, regrid_nprcy, 
                regrid_nex, regrid_ney, regrid_nez, 
                regrid_porder ): 
    conf_run_s = f"""#--- Configuration file for a test case of gravity wave  -------
&PARAM_IO
 IO_LOG_BASENAME = "regrid_LOG"
! IO_LOG_ALLNODE  = .true., 
/
&PARAM_CONST
  CONST_OHM = 0.0D0, 
/
&PARAM_REGRID_MESH
 in_MeshType  = "CUBEDSPHERE3D", 
 out_MeshType = "LONLAT3D", 
/
&PARAM_REGRID_INTERP_FIELD
  !- input --------------------
  in_basename="history",      
  vars = "W", "U", "V", "DDENS", "THERM", 
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
def get_job_header(job_name, nprc, elapse_time):
  node_num = math.ceil(nprc/4)
  if node_num > 384:
    rscgrp = "large"
  if node_num == 384:
    rscgrp = "large"
    node_num=385 
  else:
    rscgrp = "small"
  
  jobshell_s = f"""################################################################################
#
# for Fugaku
#
################################################################################
#PJM --rsc-list "rscunit=rscunit_ft01"
#PJM --name "{job_name}"
#PJM -x PJM_LLIO_GFSCACHE=/vol0005
#PJM --rsc-list "rscgrp={rscgrp}"
#PJM --rsc-list "node={node_num}"
#PJM --rsc-list "elapse={elapse_time}"
#PJM --mpi "max-proc-per-node=4"
#PJM -S


module purge
module load lang/tcsds-1.2.38

export SPACK_LIB_PATH=/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/parallel-netcdf-1.12.3-avpnzm4pwv2tuu2mv73lacb4vhcwlnds/lib:/opt/FJSVxtclanga/tcsds-mpi-latest/lib64:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/netcdf-fortran-4.6.0-mmdtg5243y4mwqsl3gcu3m2kh27raq5n/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/netcdf-c-4.9.0-g462kcd2ivou7ewax6wddywoyrbz2oib/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/hdf5-1.12.2-kb4msz2kuwzsmqsshhpryqebui6tqcfs/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/rhash-1.4.2-s3mitrsnpm36uemub4vkzj22qa4ygndu/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/libuv-1.44.1-riv7xhqvpur57jexesqfpw2mpnjjfhdd/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/libarchive-3.5.2-l7jdc7uw35jngg7tibqzsohz44ouwsj7/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/zstd-1.5.2-7j2edrlmibpft52s3m3q7ujechw3hujt/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/zlib-1.2.13-go4ye2sg72pcca4bgunmcseuzq6czbol/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/mbedtls-2.28.0-squ3v2xuqnd3mfpxiuoimtxaookk3dyi/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/lzo-2.10-uhskbd2ewdp4akltdmetra3oy4twv57f/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/libiconv-1.16-bfdxvmujixuefjz26ldcsxhzqr3rcufm/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/expat-2.4.8-lztkevt2hobbf7ykiwnuegynnoxqqvwe/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/libbsd-0.11.5-x462pikjmy4scmsuhucngco5efautbg2/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/libmd-1.0.4-wcmufmjxfiwxa65p4eetl2y674q2pgqa/lib
export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:/opt/FJSVxtclanga/tcsds-mpi-latest/lib64:${{SPACK_LIB_PATH}}:${{LD_LIBRARY_PATH}}

#export XOS_MMM_L_ARENA_FREE=1
export FORT90L="-Wl,-T"
#export OMPI_MCA_plm_ple_memory_allocation_policy=bind_local
export PLE_MPI_STD_EMPTYFILE="off"
export PARALLEL=12
export OMP_NUM_THREADS=12
#export fu11bf=1

SCALE_DG_INIT_BIN={SCALE_DG_BIN_PATH}/scale-dg_init
SCALE_DG_BIN={SCALE_DG_BIN_PATH}/scale-dg
SCALE_DG_REGRID_BIN={SCALE_DG_REGRID_BIN_PATH}/regrid_tool
  """
  return jobshell_s  

def mksh_job_run( conf_path, job_name, 
              nprc, elapse_time ):

  jobshell_header_s = get_job_header(job_name, nprc, elapse_time)
  jobshell_s = f"""
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

  jobshell_header_s = get_job_header(job_name, nprc, elapse_time)
  jobshell_s = f"""
mkdir -p {outdir}/
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
    porder = exp_info["porder"]

    out_dir_pref=f"./rhot_hevi/{exp_name}"

    print(out_dir_pref)
    os.makedirs(out_dir_pref, exist_ok=True)
    
    mkconf_init(f"{out_dir_pref}/init.conf", 
                nprc, eh, ez, porder, exp_info["initgp_porder"] )
    
    mkconf_run(f"{out_dir_pref}/run.conf", 
               "init_00000101-000000.000", 1, 
                nprc, eh, ez, porder, 
                exp_info["dt"], exp_info["dt_dyn"],
                exp_info["modal_filter_flag"], exp_info["mf_alph"], exp_info["mf_ordh"], exp_info["mf_alpv"], exp_info["mf_ordv"])        
                    
    mkconf_regrid(f"{out_dir_pref}/regrid.conf", 
                    nprc, eh, ez, porder,  
                    exp_info["regrid_nprcx"], exp_info["regrid_nprcy"], 
                    exp_info["regrid_Ex"], exp_info["regrid_Ey"], exp_info["Ez"], exp_info["regrid_porder"] )
                
    mksh_job_run(f"{out_dir_pref}/job_run.sh", f"GRAV_Eh{eh}P{porder}", 
                      nprc, exp_info["elapse_time"]) 
        
    mksh_job_regrid(f"{out_dir_pref}/job_regrid.sh", f"REG_E{eh}P{porder}", "regrid.conf", 
                      exp_info["regrid_nprcx"]*exp_info["regrid_nprcy"], exp_info["regrid_elapse_time"], 
                      "outdata")
            
    mksh_job_regrid(f"{out_dir_pref}/job_regrid_topo.sh", f"REGT_E{eh}P{porder}", "regrid_topo.conf", 
                      exp_info["regrid_nprcx"]*exp_info["regrid_nprcy"], exp_info["regrid_elapse_time"], 
                      "outdata")
  
#---------------------------------
