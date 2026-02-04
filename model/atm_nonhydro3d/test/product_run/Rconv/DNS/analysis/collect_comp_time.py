import os
import sys
import numpy as np
sys.path.append(os.path.join(os.path.dirname('__file__'), '../../analysis_lib/'))
import comp_time_analysis_common as common

#--
OUT_DIR = "./comp_cost_info/"
EXP_DIR="../"

exp_list = {
  "Dx25m_P3": { "run_s": 4, "run_e": 20 }, 
  "Dx25m_P7": { "run_s": 4, "run_e": 20 }, 
  "Dx27m_P11": { "run_s": 4, "run_e": 20 }, 
  "Dx12.5m_P3": { "run_s": 8, "run_e": 24 }, 
  "Dx12.5m_P7": { "run_s": 8, "run_e": 24 }, 
  "Dx13m_P11": { "run_s": 8, "run_e": 24 },   
  "Dx6.3m_P3": { "run_s": 13, "run_e": 31 },   
  "Dx6.3m_P7": { "run_s": 8, "run_e": 26 },   
  "Dx6.7m_P11": { "run_s": 13, "run_e": 34 },     
}

node_info_list = {
    "Dx25m_P3": 16, "Dx25m_P7": 16, "Dx27m_P11": 25, 
    "Dx12.5m_P3": 64, "Dx12.5m_P7": 64, "Dx13m_P11": 100,
    "Dx6.3m_P3": 256*2, "Dx6.3m_P7": 256*2, "Dx6.7m_P11":400*1,  # x2 for P3, P7 becaues integration time is 3600 s per run
    "Dx6.3m_P3_comp1": 256*4, "Dx6.3m_P7_comp1": 256*4, "Dx6.7m_P11_comp1":400*4, # x4 becaues integration time is 1800 s per run
}

sect_list = [
    'MAIN ATM_DYN_core', 
    'MAIN ATM_DYN_update_caltend_ex', 
    # 'MAIN ATM_DYN_exchange_prgv',
    'MAIN ATM_DYN_exchange_prgv_wait',    
    'MAIN ATM_DYN_update_advance',
    'MAIN ATM_DYN_applyBC_prgv', 
    'MAIN ATM_DYN_cal_pres',
    'MAIN ATM_DYN_update_modalfilter',
    'MAIN ATM_DYN_update_add_tp', 
    'MAIN ATM_DYN_update_pre',
    'MAIN ATM_DYN_update_post'    
]

#-
def get_prof_rapreport( exp_name, run_s, run_e ):
  prof_info = []
  for runno in range(run_s,run_e+1):
      prof_info_ = common.get_prof_rapreport(f"{EXP_DIR}/{exp_name}/run{runno}/LOG.pe000000")
      prof_info.append(prof_info_)
  return prof_info

#--
prof_info_list = {}
for exp_name, exp_info in exp_list.items():
    prof_info_list[exp_name] = get_prof_rapreport(exp_name, exp_info["run_s"], exp_info["run_e"])


if not os.path.exists(OUT_DIR):
   os.makedirs(OUT_DIR)

for exp_name in exp_list.keys():
  print(f"* {exp_name} ------------------------------------")
  common.output_comp_time(prof_info_list[exp_name], f"{OUT_DIR}/{exp_name}.dat")

common.output_comp_resource_mean_usage(prof_info_list, sect_list, node_info_list, True, 
                                       f"{OUT_DIR}/comp_resource_usage.dat" )