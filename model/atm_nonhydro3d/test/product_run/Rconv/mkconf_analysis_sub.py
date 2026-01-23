import sys
import os
import math
import datetime
sys.path.append(os.path.join(os.path.dirname(__file__), './common'))
import batch_job_common

SCALE_DG_BIN_PATH="../../"
SCALE_DG_REGRID_BIN_PATH="../../../../../../../../bin"

#----------------------

def mkconf_analysis( conf_path,
                nprcx, nex, nprcy, ney, nez, porder, 
                in_bs_filebase, 
                start_time_sec, out_tintrv, num_step, tb_params, 
                bc_fixed_heat_flux, bc_stabcoef ):
    
    tb_scheme = tb_params["type"]
    if tb_scheme == "DNS":
        dns_mu = tb_params["dns_mu"]; dns_nu = tb_params["dns_nu"]; 
        tb_params = f"""
&PARAM_ATMOS_PHY_TB_DGM_DNS
  DNS_MU       = {dns_mu}, 
  DNS_NU       = {dns_nu}, 
/
"""
    else:
      tb_scheme = "SMAGORINSKY" 
      tb_params = ""
    
    conf_init_s = f"""#--- Configuration file for analysis  -------
&PARAM_IO
 IO_LOG_BASENAME = "analysis_LOG"
! IO_LOG_ALLNODE  = .true., 
/
&PARAM_RBCONV_ANALYSIS
  in_filebase      = "history",
  in_bs_filebase   = "{in_bs_filebase}",
  out_filebase_V1D = "analysis/v1D", 
  out_filebase_tb = "analysis/tb",
  start_time0    = {start_time_sec}, 
  output_tintrv  = {out_tintrv},
  num_step       = {num_step}, 
  tb_scheme      = "{tb_scheme}", 
  U0             = 5D0, 
/
&PARAM_RBCONV_ANALYSIS_BC
  BTM_BC_TYPE_HEAT = 'FixedFlux', 
  TOP_BC_TYPE_HEAT = 'FixedFlux', 
  BTM_FIXED_HEAT_FLUX = {bc_fixed_heat_flux}, 
  TOP_FIXED_HEAT_FLUX = {bc_fixed_heat_flux}, 
  StabCoef_bnd        = {bc_stabcoef}, 
/
&PARAM_ATMOS_MESH
  dom_xmin         =   -1.6D3,  
  dom_xmax         =    1.6D3, 
  dom_ymin         =   -1.6D3,  
  dom_ymax         =    1.6D3,  
  dom_zmin         =   0.0D0,  
  dom_zmax         =   1.6D3, 
  NprcX            = {nprcx}, 
  NeX              = {nex},
  NprcY            = {nprcy}, 
  NeY              = {ney},  
  NeZ              = {nez}, 
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},
/ 
&PARAM_ATMOS_DYN_BND
  btm_vel_bc   = 'NOSLIP', 
  top_vel_bc   = 'NOSLIP', 
  btm_thermal_bc   = 'ADIABATIC', 
  top_thermal_bc   = 'ADIABATIC', 
  north_thermal_bc  = 'ADIABATIC',
  south_thermal_bc  = 'ADIABATIC', 
  west_thermal_bc   = 'ADIABATIC',
  east_thermal_bc   = 'ADIABATIC',
/
{tb_params}
    """
    
    with open(conf_path, 'w') as f:
        f.write(conf_init_s)

#----------------
def mksh_job_analysis( conf_path, job_name, conf_name, 
              nprc, elapse_time, outdir ):

  jobshell_header_s = batch_job_common.get_job_header(job_name, nprc, elapse_time)
  jobshell_s = f"""
SCALE_DG_ANALYSIS_BIN={SCALE_DG_REGRID_BIN_PATH}/RBconv_analysis_fr
mkdir -p {outdir}/  
llio_transfer ${{SCALE_DG_ANALYSIS_BIN}} *.conf

mpiexec -np {nprc} -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr \\
  ${{SCALE_DG_ANALYSIS_BIN}} {conf_name} || exit 1

llio_transfer --purge ${{SCALE_DG_ANALYSIS_BIN}} *.conf
  """
  
  with open(conf_path, 'w') as f:
      f.write(jobshell_header_s + jobshell_s)


#----------------
def mk_conf_sh( exp_name, exp_info, tb_params ):
    porder = exp_info["porder"]         
    exp_top_dir = exp_info["top_dir"]
    analysis_run_no_s = exp_info["analysis_run_no_s"]
    analysis_run_no_e = exp_info["analysis_run_no_e"]
    analysis_start_time_sec = exp_info["analysis_start_time_sec"]
    analysis_integ_time_per_run = exp_info["analysis_integ_time_per_run"]
    analysis_out_tintrv = exp_info["analysis_out_tintrv"]
    
    
    for runno in range(analysis_run_no_s, analysis_run_no_e+1):
      
      out_dir_pref=f"./{exp_top_dir}/{exp_name}/run{runno}"
      print(out_dir_pref)
      # os.makedirs(out_dir_pref, exist_ok=True)
      start_time_sec = analysis_start_time_sec + (runno - analysis_run_no_s) * analysis_integ_time_per_run
      num_step = int(analysis_integ_time_per_run / analysis_out_tintrv) + 1
      
      
      mkconf_analysis(f"{out_dir_pref}/analysis.conf", 
                      exp_info["nprcx"], exp_info["Ex"], exp_info["nprcy"], exp_info["Ey"], exp_info["Ez"], porder, 
                      exp_info["analysis_in_bs_filebase"], start_time_sec, analysis_out_tintrv, num_step, tb_params,
                      exp_info["const_hflx"], exp_info["StabCoef_bnd"] )
                  
      mksh_job_analysis(f"{out_dir_pref}/job_analysis.sh", f"RB_ANALYSIS", "analysis.conf", 
                        exp_info["nprcx"]*exp_info["nprcy"], exp_info["analysis_elapse_time"], 
                        "analysis")
  
#---------------------------------
