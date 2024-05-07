import os
import math
import datetime
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../common'))
import batch_job_common


SCALE_DG_BIN_PATH="../../../"
SCALE_DG_REGRID_BIN_PATH="../../../../../../../../bin"

#----------------------
def mkconf_sh_spectra( conf_path,
                nprc, neh, nez, porder, fz, rplanet, shallow_atm_approx, 
                Mt, out_dir ):

  out_dir="./analysis"    
     
  if shallow_atm_approx:
    shallow_atm_approx_flag = ""
  else:
    shallow_atm_approx_flag = "SHALLOW_ATM_APPROX_FLAG = .false.,"
     
  conf_run_s = f"""#--- Configuration file for SH transform tools  -------
&PARAM_IO
 IO_LOG_BASENAME = "SH_LOG"
! IO_LOG_ALLNODE  = .true., 
/  
&PARAM_CONST
  CONST_RADIUS = {rplanet}D0,  
/
&PARAM_SH_TRANSFORM
  in_filebase="history",  
  out_filebase="{out_dir}/spectral_data",  
  Mt                  = {Mt},
  LevelNum            = 1, 
  TARGET_LEVELS       = 500D0, 
  LEVEL_UNITS         = "m",
  VARS                = "U", "V", "W", "PT", 
  TARGET_PROC_NUM_TOT = {nprc}, 
  KinEnergyAnalysisFlag = .true., 
/
&PARAM_SH_MESH
  !---------------------
  {shallow_atm_approx_flag}  
  Nprc             = {nprc}, 
  NeGX             = {neh},
  NeGY             = {neh},
  NeZ              = {nez},
  dom_zmin         = 0.0D0, 
  dom_zmax         = 3.0D3,  
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},
  FZ               = {fz},  
/
    """
    
  with open(conf_path, 'w') as f:
      f.write(conf_run_s)

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
    hist_int_sec = exp_info["hist_int_sec_analysis"] 
    
    for runno in range(1,run_num+1):
      date_time = date_time0 + datetime.timedelta(hours=hr_per_run*(runno-1))
      
      out_dir_pref = f"{out_dir_pref0}/run{runno}/"
      os.makedirs(out_dir_pref, exist_ok=True)
      

      mkconf_sh_spectra(f"{out_dir_pref}/spectral_analysis.conf", 
                      nprc, eh, ez, porder, fz, rplanet, shallow_atm_approx,
                      exp_info["spectra_analysis_Mt"], "analysis" )
                  
      mksh_job_spectra_analysis(f"{out_dir_pref}/job_spectra_analysis.sh", f"SA_E{eh}P{porder}", "spectral_analysis.conf", 
                        exp_info["spectra_analysis_nprc"], exp_info["spectra_analysis_elapse_time"], 
                        "analysis")
  
#---------------------------------
