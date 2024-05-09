import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../common'))
import batch_job_common

OUTPUT_INT=600
OUTPUT_NSTEP=13

SCALE_DG_MTWAVE_ANALYSIS_BIN="../../../../../../../bin/mtwave_analysis_rm"
SCALE_DG_REGRID_BIN_PATH="../../../../../../../bin"

#----------------------

def mkconf_run_analysis( conf_path, 
                nprc, neh, nez, porder, fz, 
                in_ref_fbase, in_bs_ref_fbase, in_topo_ref_fbase, 
                dom_xmin, dom_xmax, dom_dy, 
                window_xs, window_xe, window_zs, window_ze, 
                start_time, dt, nstep, outdata_analysis, outdata_bs_analysis, analysis_dir ): 
    conf_run_s = f"""#--- Configuration file for analysis  -------
&PARAM_IO
 IO_LOG_BASENAME = 'analysis_LOG',
/
&PARAM_MTWAVE_ANALYSIS
  in_filebase ='{outdata_analysis}/history',
  in_bs_filebase ='{outdata_bs_analysis}/bs',
  in_ref_filebase ='{in_ref_fbase}',
  in_bs_ref_filebase ='{in_bs_ref_fbase}', 
  in_topo_ref_filebase ='{in_topo_ref_fbase}', 
!  out_filebase_V1D='analysis/v1D',
  start_time0   = {start_time}D0,
  output_tintrv = {dt}D0,
  num_step            = {nstep}, 
  PolyOrderErrorCheck = 11, 
  NUMERROR_LOG_OUT_BASENAME = '{analysis_dir}/NUMERROR_LOG', 
  window_xs = {window_xs}, 
  window_xe = {window_xe}, 
  window_zs = {window_zs}, 
  window_ze = {window_ze}, 
/
&PARAM_ATMOS_MESH
  NprcX            = {nprc}, 
  NprcY            = 1,   
  NeX              = {neh},
  NeY              = 1,
  NeZ              = {nez}, 
  dom_xmin         = {dom_xmin}, 
  dom_xmax         = {dom_xmax},
  dom_ymin         = 0.0D3,  
  dom_ymax         = {dom_dy}, 
  dom_zmin         = 0.0D0, 
  dom_zmax         = 30.0D3, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},
  FZ               = {fz}, 
  LumpedMassMatFlag = .false.,
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

SCALE_DG_MTWV_ANALYSIS_BIN={SCALE_DG_MTWAVE_ANALYSIS_BIN}
llio_transfer ${{SCALE_DG_MTWV_ANALYSIS_BIN}} *.conf

mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_MTWV_ANALYSIS_BIN}} analysis.conf || exit 1  

llio_transfer --purge ${{SCALE_DG_MTWV_ANALYSIS_BIN}} *.conf    
  """
  
  with open(conf_path, 'w') as f:
      f.write(jobshell_header_s + jobshell_s)

             
def mk_conf_sh( exp_name, exp_info, out_dir_pref0="./rhot_heve_v2" ):
    nprcx = exp_info["nprc"]
    eh = exp_info["Eh"]
    dom_xmin = exp_info["dom_xmin"]
    dom_xmax = exp_info["dom_xmax"]
    dom_dy = exp_info["dom_dy"]        
    porder = exp_info["porder"]
    regrid_Nprc = exp_info["regrid_nprc"]    
    regrid_Eh = exp_info["regrid_Ex"]
    regrid_Ez = exp_info["regrid_Ez"]
    regrid_Porder = exp_info["regrid_porder"]
    regrid_FZ = exp_info["regrid_fz"]
    refsol_dir = exp_info["refsol_dir"]
    

    out_dir=f"{out_dir_pref0}/{exp_name}"

    print(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    
    mkconf_run_analysis(f"{out_dir}/analysis.conf", 
                regrid_Nprc, regrid_Eh, regrid_Ez, regrid_Porder, regrid_FZ, 
                f'{refsol_dir}/outdata_compari/history',
                f'{refsol_dir}/outdata_compari/init_00000101-000000.000',
                f'{refsol_dir}/outdata_compari/TOPO', 
                dom_xmin, dom_xmax, dom_dy, 
                exp_info["window_xs"], exp_info["window_xe"], exp_info["window_zs"], exp_info["window_ze"],   
                0, OUTPUT_INT, OUTPUT_NSTEP, 
                "outdata_compari", "outdata_compari", "analysis")        
                                     
    mksh_job_analysis( f"{out_dir}/job_analysis.sh", f"ANL_E{nprcx*eh}P{porder}", 
                regrid_Nprc, exp_info["regrid_elapse_time"], "analysis" )
    
  
#---------------------------------