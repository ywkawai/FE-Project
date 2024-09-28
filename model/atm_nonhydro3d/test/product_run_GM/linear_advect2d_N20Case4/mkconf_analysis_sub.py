import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../common'))
import batch_job_common

REGRID_Nprc=1536 
REGRID_Eh=64
REGRID_Ez=3
REGRID_Porder=7
REFSOL_DIR_SUFFIX="_check"

OUTPUT_INT=86400
OUTPUT_NSTEP=13

SCALE_DG_BIN_PATH="../.."
SCALE_DG_REGRID_BIN_PATH="../../../../../../../bin"
SCALE_DG_LINADV_ANALYSIS_BIN="../../../../../../../bin/linadv2d_N20Case4_analysis"

#----------------------

def mkconf_run_analysis( conf_path, 
                nprc, neh, nez, porder, 
                in_ref_fbase, 
                start_time, dt, nstep, outdata_analysis, analysis_dir ): 
    conf_run_s = f"""#--- Configuration file for analysis  -------
&PARAM_IO
 IO_LOG_BASENAME = 'analysis_LOG',
/
&PARAM_LINADV_ANALYSIS
  in_filebase ='{outdata_analysis}/history',
  in_ref_filebase ='{in_ref_fbase}',
  start_time0   = {start_time}D0,
  output_tintrv = {dt}D0,
  num_step            = {nstep}, 
  PolyOrderErrorCheck = 11, 
  NUMERROR_LOG_OUT_BASENAME = '{analysis_dir}/NUMERROR_LOG', 
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
    conf_run_s = f"""#--- Configuration file for a test case of mountain wave  -------
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

def mksh_job_analysis( conf_path, job_name, 
              nprc, elapse_time, outdir ):

  jobshell_header_s = batch_job_common.get_job_header(job_name, nprc, elapse_time)
  jobshell_s = f"""
mkdir -p {outdir}/

SCALE_DG_LINADV_ANALYSIS_BIN={SCALE_DG_LINADV_ANALYSIS_BIN}
llio_transfer ${{SCALE_DG_LINADV_ANALYSIS_BIN}} *.conf

mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_LINADV_ANALYSIS_BIN}} analysis.conf || exit 1  

llio_transfer --purge ${{SCALE_DG_LINADV_ANALYSIS_BIN}} *.conf    
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
            
  jobshell_header_s = batch_job_common.get_job_header(job_name, nprc, elapse_time)
  jobshell_s = f"""
SCALE_DG_REGRID_BIN={SCALE_DG_REGRID_BIN_PATH}/regrid_tool

mkdir -p {outdir}/
llio_transfer ${{SCALE_DG_REGRID_BIN}} *.conf

{run_cmd}

llio_transfer --purge ${{SCALE_DG_REGRID_BIN}} *.conf 
  """
  
  with open(conf_path, 'w') as f:
      f.write(jobshell_header_s + jobshell_s)
             
def mk_conf_sh( ptracer_shape, exp_name, exp_info ):
    nprc = exp_info["nprc"]
    eh = exp_info["Eh"]
    ez = exp_info["Ez"]
    porder = exp_info["porder"]

    out_dir=f"./{ptracer_shape}/{exp_name}"

    print(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    
    refsol_dir=f"../Eh{REGRID_Eh}Ez{REGRID_Ez}P{REGRID_Porder}{REFSOL_DIR_SUFFIX}"
    mkconf_run_analysis(f"{out_dir}/analysis.conf", 
                REGRID_Nprc, REGRID_Eh, REGRID_Ez, REGRID_Porder, 
                f'{refsol_dir}/history',
                0, OUTPUT_INT, OUTPUT_NSTEP, 
                "outdata_analysis", "analysis")        
                                     
    mksh_job_analysis( f"{out_dir}/job_analysis.sh", f"ANL_E{eh}P{porder}", 
                REGRID_Nprc, exp_info["regrid_elapse_time"], "analysis" )
    
    mkconf_regrid(f"{out_dir}/regrid_analysis.conf", 
                nprc, eh, ez, porder, 
                REGRID_Nprc, REGRID_Eh, REGRID_Ez, REGRID_Porder, 
                "regrid_analysis_LOG", "history", "'PTracer'", "./outdata_analysis/history" )
    
    regrid_conf_list = [ "regrid_analysis.conf" ]
    regrid_conf_list.append("regrid_bs_analysis.conf")
    mksh_job_regrid( f"{out_dir}/job_regrid_analysis.sh", f"REG_E{eh}P{porder}", 
                regrid_conf_list, REGRID_Nprc, exp_info["regrid_elapse_time"], "outdata_analysis")
  
#---------------------------------