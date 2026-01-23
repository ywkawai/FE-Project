import mkconf_analysis_sub as mkconf_sub

exp_info_common = {"top_dir": "LES",
                   "const_hflx": "200D0", "StabCoef_bnd": "1.0D0" }
tb_params = {"type": "LES", "dns_nu": "0D0", "dns_mu": "0D0",}

# p=3
exp_list_P3 = {
  "Dx25m_P3": { "nprcx": 8, "nprcy": 8, "Ex": 4, "Ey": 4, "Ez": 16, "porder": 3, "analysis_elapse_time": "00:10:00", 
                  "analysis_run_no_s": 3, "analysis_run_no_e": 24, 
                  "analysis_start_time_sec": 2*7200, "analysis_integ_time_per_run": 7200.0,  "analysis_out_tintrv": 900.0,
                  "analysis_in_bs_filebase": "../run2/restart_00010101-040000.000", 
                },    
  "Dx12.5m_P3": { "nprcx": 16, "nprcy": 16,  "Ex": 4, "Ey": 4, "Ez": 32, "porder": 3, "analysis_elapse_time": "00:10:00", 
                  "analysis_run_no_s": 3, "analysis_run_no_e": 24, 
                  "analysis_start_time_sec": 2*7200, "analysis_integ_time_per_run": 7200.0,  "analysis_out_tintrv": 900.0,
                  "analysis_in_bs_filebase": "../run2/restart_00010101-040000.000", 
                },          
}

# p=7
exp_list_P7 = {
  "Dx25m_P7": { "nprcx": 8, "nprcy": 8, "Ex": 2, "Ey": 2, "Ez": 8, "porder": 7,  "analysis_elapse_time": "00:5:00", 
                  "analysis_run_no_s": 3, "analysis_run_no_e": 24, 
                  "analysis_start_time_sec": 2*7200, "analysis_integ_time_per_run": 7200.0,  "analysis_out_tintrv": 900.0,
                  "analysis_in_bs_filebase": "../run1/init_00010101-000000.000", 
                },    
  "Dx12.5m_P7": { "nprcx": 16, "nprcy": 16, "Ex": 2, "Ey": 2, "Ez": 16, "porder": 7,  "analysis_elapse_time": "00:10:00", 
                  "analysis_run_no_s": 3, "analysis_run_no_e": 24, 
                  "analysis_start_time_sec": 2*7200, "analysis_integ_time_per_run": 7200.0,  "analysis_out_tintrv": 900.0,
                  "analysis_in_bs_filebase": "../run1/init_00010101-000000.000", 
                },
}
# p=11
exp_list_P11 = {
  "Dx27m_P11": { "nprcx": 10, "nprcy": 10, "Ex": 1, "Ey": 1, "Ez": 5, "porder": 11,  "analysis_elapse_time": "00:10:00", 
                  "analysis_run_no_s": 3, "analysis_run_no_e": 24, 
                  "analysis_start_time_sec": 2*7200, "analysis_integ_time_per_run": 7200.0,  "analysis_out_tintrv": 900.0,
                  "analysis_in_bs_filebase": "../run2/restart_00010101-040000.000", 
                },   
  "Dx13m_P11": { "nprcx": 20, "nprcy": 20, "Ex": 1, "Ey": 1, "Ez": 10, "porder": 11,  "analysis_elapse_time": "00:10:00", 
                  "analysis_run_no_s": 3, "analysis_run_no_e": 24, 
                  "analysis_start_time_sec": 2*7200, "analysis_integ_time_per_run": 7200.0,  "analysis_out_tintrv": 900.0,
                  "analysis_in_bs_filebase": "../run2/restart_00010101-040000.000", 
                }, 
}

exp_list = {  
           **exp_list_P3,                          
           **exp_list_P7,  
           **exp_list_P11,              
            }
#---------------------------------

for exp_name, exp_info in exp_list.items():
  exp_info.update(exp_info_common)
  mkconf_sub.mk_conf_sh( exp_name, exp_info, tb_params )
