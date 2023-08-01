import os
import math

REGRID_Nprc=1536
REGRID_Eh=64
REGRID_Ez=48
REGRID_Porder=7

OUTPUT_INT=4*3600
OUTPUT_NSTEP=13

SCALE_DG_GRAVWAVE_ANALYSIS_BIN="../../../../../../../bin/gravwave_analysis"
SCALE_DG_REGRID_BIN_PATH="../../../../../../../bin"

#----------------------

def mkconf_run_analysis( conf_path, 
                nprc, neh, nez, porder, 
                in_ref_fbase, in_bs_ref_fbase,  
                start_time, dt, nstep, outdata_analysis, outdata_bs_analysis, analysis_dir ): 
    conf_run_s = f"""#--- Configuration file for analysis  -------
&PARAM_IO
 IO_LOG_BASENAME = 'analysis_LOG',
/
&PARAM_GRAVWAVE_ANALYSIS
  in_filebase='{outdata_analysis}/history',
  in_bs_filebase='{outdata_bs_analysis}/bs',
  in_ref_filebase='{in_ref_fbase}',
  in_bs_ref_filebase='{in_bs_ref_fbase}', 
!  out_filebase_V1D='analysis/v1D',
  start_time0={start_time}D0,
  output_tintrv={dt}D0,
  num_step={nstep}, 
  PolyOrderErrorCheck=11, 
  NUMERROR_LOG_OUT_BASENAME='{analysis_dir}/NUMERROR_LOG', 
/
&PARAM_ATMOS_MESH
  Nprc             = {nprc}, 
  NeGX             = {neh},
  NeGY             = {neh},
  NeZ              = {nez},
  dom_zmin         = 0.0D0, 
  dom_zmax         = 10.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},
  LumpedMassMatFlag = .false.,
/
    """
    
    with open(conf_path, 'w') as f:
        f.write(conf_run_s)

def mkconf_regrid( conf_path,
                nprc, neh, nez, porder, 
                regrid_nprc, regrid_neh, regrid_nez, regrid_porder, 
                log_basename, in_basename, vars_list, regrid_outdata ): 
    conf_run_s = f"""#--- Configuration file for a test case of sound wave  -------
&PARAM_IO
 IO_LOG_BASENAME = "{log_basename}", 
! IO_LOG_ALLNODE  = .true., 
/
&PARAM_REGRID_MESH
 in_MeshType  = "CUBEDSPHERE3D", 
 out_MeshType = "CUBEDSPHERE3D", 
/
&PARAM_REGRID_INTERP_FIELD
  !- input --------------------
  in_basename="{in_basename}",      
  vars = {vars_list}, 
  !out_tinterval = 5,
/
&PARAM_REGRID_FILE
  !-- output ----------------
  out_basename="{regrid_outdata}", 
  out_UniformGrid=.false., 
  out_dtype="REAL8", 
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
&PARAM_REGRID_OUTMESH3D_CUBEDSPHERE
  NLocalMeshPerPrc = 1, 
  Nprc             = {regrid_nprc}, 
  NeGX             = {regrid_neh},
  NeGY             = {regrid_neh},
  NeGZ             = {regrid_nez},
  dom_zmin         = 0.0D0, 
  dom_zmax         = 10.0D3, 
  PolyOrder_h      = {regrid_porder},
  PolyOrder_v      = {regrid_porder},
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
    node_num = 385
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
module load lang/tcsds-1.2.37

export SPACK_LIB_PATH=/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/parallel-netcdf-1.12.3-avpnzm4pwv2tuu2mv73lacb4vhcwlnds/lib:/opt/FJSVxtclanga/tcsds-mpi-latest/lib64:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/netcdf-fortran-4.6.0-mmdtg5243y4mwqsl3gcu3m2kh27raq5n/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/netcdf-c-4.9.0-g462kcd2ivou7ewax6wddywoyrbz2oib/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/hdf5-1.12.2-kb4msz2kuwzsmqsshhpryqebui6tqcfs/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/rhash-1.4.2-s3mitrsnpm36uemub4vkzj22qa4ygndu/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/libuv-1.44.1-riv7xhqvpur57jexesqfpw2mpnjjfhdd/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/libarchive-3.5.2-l7jdc7uw35jngg7tibqzsohz44ouwsj7/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/zstd-1.5.2-7j2edrlmibpft52s3m3q7ujechw3hujt/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/zlib-1.2.13-go4ye2sg72pcca4bgunmcseuzq6czbol/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/mbedtls-2.28.0-squ3v2xuqnd3mfpxiuoimtxaookk3dyi/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/lzo-2.10-uhskbd2ewdp4akltdmetra3oy4twv57f/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/libiconv-1.16-bfdxvmujixuefjz26ldcsxhzqr3rcufm/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/expat-2.4.8-lztkevt2hobbf7ykiwnuegynnoxqqvwe/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/libbsd-0.11.5-x462pikjmy4scmsuhucngco5efautbg2/lib:/vol0004/apps/oss/spack-v0.19/opt/spack/linux-rhel8-a64fx/fj-4.8.1/libmd-1.0.4-wcmufmjxfiwxa65p4eetl2y674q2pgqa/lib
export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:/opt/FJSVxtclanga/tcsds-mpi-latest/lib64:${{SPACK_LIB_PATH}}:${{LD_LIBRARY_PATH}}

#export XOS_MMM_L_ARENA_FREE=1
export FORT90L="-Wl,-T"
#export OMPI_MCA_plm_ple_memory_allocation_policy=bind_local
export PLE_MPI_STD_EMPTYFILE="off"
export PARALLEL=12
export OMP_NUM_THREADS=12
#export fu11bf=1

SCALE_DG_GRAVWAVE_ANALYSIS_BIN={SCALE_DG_GRAVWAVE_ANALYSIS_BIN}
SCALE_DG_REGRID_BIN={SCALE_DG_REGRID_BIN_PATH}/regrid_tool
  """
  return jobshell_s  

def mksh_job_analysis( conf_path, job_name, 
              nprc, elapse_time, outdir ):

  jobshell_header_s = get_job_header(job_name, nprc, elapse_time)
  jobshell_s = f"""
mkdir -p {outdir}/
llio_transfer ${{SCALE_DG_GRAVWAVE_ANALYSIS_BIN}} *.conf

mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_GRAVWAVE_ANALYSIS_BIN}} analysis.conf || exit 1  

llio_transfer --purge ${{SCALE_DG_GRAVWAVE_ANALYSIS_BIN}} *.conf    
  """
  
  with open(conf_path, 'w') as f:
      f.write(jobshell_header_s + jobshell_s)

#----------------
def mksh_job_regrid( conf_path, job_name, regrid_cnf_list, 
              nprc, elapse_time, outdir ):

  run_cmd_list = []
  for regrid_cnf in regrid_cnf_list:
    run_cmd_list.append(f"""mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
      ${{SCALE_DG_REGRID_BIN}} {regrid_cnf} || exit 1
                   """)
    
  run_cmd = "\n".join(run_cmd_list)
            
  jobshell_header_s = get_job_header(job_name, nprc, elapse_time)
  jobshell_s = f"""
mkdir -p {outdir}/
llio_transfer ${{SCALE_DG_REGRID_BIN}} *.conf

{run_cmd}

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
        
    refsol_dir_pref=f"../Eh{REGRID_Eh}Ez{REGRID_Ez}P{REGRID_Porder}"
    mkconf_run_analysis(f"{out_dir_pref}/analysis.conf", 
                REGRID_Nprc, REGRID_Eh, REGRID_Ez, REGRID_Porder, 
                f'{refsol_dir_pref}_1/history',
                f'{refsol_dir_pref}_1/init_00000101-000000.000',
                0, OUTPUT_INT, OUTPUT_NSTEP, 
                "outdata_analysis", "outdata_analysis", "analysis")        
                      
    mkconf_regrid(f"{out_dir_pref}/regrid_analysis.conf", 
                nprc, eh, ez, porder, 
                REGRID_Nprc, REGRID_Eh, REGRID_Ez, REGRID_Porder, 
                "regrid_analysis_LOG", "history", """THERM""", "./outdata_analysis/history" )

    mkconf_regrid(f"{out_dir_pref}/regrid_bs_analysis.conf", 
                nprc, eh, ez, porder, 
                REGRID_Nprc, REGRID_Eh, REGRID_Ez, REGRID_Porder, 
                "regrid_bs_analysis_LOG", "init_00000101-000000.000", "'PRES_hyd', 'DENS_hyd'", "./outdata_analysis/bs" )
    
    regrid_conf_list = [ "regrid_analysis.conf", "regrid_bs_analysis.conf" ]
    mksh_job_regrid(f"{out_dir_pref}/job_regrid_analysis.sh", f"REG_E{eh}P{porder}", 
                regrid_conf_list, REGRID_Nprc, exp_info["regrid_elapse_time"], "outdata_analysis")
  
#---------------------------------