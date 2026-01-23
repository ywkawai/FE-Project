import mkconf_analysis_sub as mkconf_sub

exp_info_common = {"top_dir": "DNS",
                   "const_hflx": "200D0", "StabCoef_bnd": "7.5D0" }
tb_params = {"type": "DNS", "dns_nu": "4D-1", "dns_mu": "5.714D-1"}

# p=3
exp_list_P3 = {
  "Dx25m_P3": { "nprcx": 8, "nprcy": 8, "Ex": 4, "Ey": 4, "Ez": 16, "porder": 3, "analysis_elapse_time": "00:10:00", 
                  "analysis_run_no_s": 3, "analysis_run_no_e": 24, 
                  "analysis_start_time_sec": 14400, "analysis_integ_time_per_run": 7200.0,  "analysis_out_tintrv": 900.0,
                  "analysis_in_bs_filebase": "../run2/restart_00010101-040000.000", 
                },    
  "Dx12.5m_P3": { "nprcx": 16, "nprcy": 16,  "Ex": 4, "Ey": 4, "Ez": 32, "porder": 3, "analysis_elapse_time": "00:10:00", 
                  "analysis_run_no_s": 8, "analysis_run_no_e": 24, 
                  "analysis_start_time_sec": 7*7200, "analysis_integ_time_per_run": 7200.0,  "analysis_out_tintrv": 900.0,
                  "analysis_in_bs_filebase": "../run8/restart_00010101-160000.000", 
                },          
  "Dx6.3m_P3": { "nprcx": 32, "nprcy": 32,  "Ex": 4, "Ey": 4, "Ez": 64, "porder": 3, "analysis_elapse_time": "00:10:00", 
                  "analysis_run_no_s": 13, "analysis_run_no_e": 57, 
                  "analysis_start_time_sec": 19*3600, "analysis_integ_time_per_run": 3600.0,  "analysis_out_tintrv": 900.0,
                  "analysis_in_bs_filebase": "../run12/restart_00010101-190000.000", 
                },              
}

# p=7
exp_list_P7 = {
  "Dx25m_P7": { "nprcx": 8, "nprcy": 8, "Ex": 2, "Ey": 2, "Ez": 8, "porder": 7,  "analysis_elapse_time": "00:5:00", 
                  "analysis_run_no_s": 3, "analysis_run_no_e": 24, 
                  "analysis_start_time_sec": 7200*2, "analysis_integ_time_per_run": 7200.0,  "analysis_out_tintrv": 900.0,
                  "analysis_in_bs_filebase": "../run1/init_00010101-000000.000", 
                },    
  "Dx12.5m_P7": { "nprcx": 16, "nprcy": 16, "Ex": 2, "Ey": 2, "Ez": 16, "porder": 7,  "analysis_elapse_time": "00:10:00", 
                  "analysis_run_no_s": 3, "analysis_run_no_e": 24, 
                  "analysis_start_time_sec": 14400, "analysis_integ_time_per_run": 7200.0,  "analysis_out_tintrv": 900.0,
                  "analysis_in_bs_filebase": "../run1/init_00010101-000000.000", 
                },
  "Dx6.3m_P7": { "nprcx": 32, "nprcy": 32,  "Ex": 2, "Ey": 2, "Ez": 32, "porder": 7, "analysis_elapse_time": "00:10:00", 
                  "analysis_run_no_s": 8, "analysis_run_no_e": 52, 
                  "analysis_start_time_sec": 7*7200, "analysis_integ_time_per_run": 3600.0,  "analysis_out_tintrv": 900.0,
                  "analysis_in_bs_filebase": "../run7/restart_00010101-140000.000", 
                },                
  "Dx3.1m_P7": { "nprcx": 64, "nprcy": 64,  "Ex": 2, "Ey": 2, "Ez": 64, "porder": 7, "analysis_elapse_time": "00:10:00", 
                  "analysis_run_no_s": 13, "analysis_run_no_e": 100, 
                  "analysis_start_time_sec": 19*3600, "analysis_integ_time_per_run": 1800.0,  "analysis_out_tintrv": 900.0,
                  "analysis_in_bs_filebase": "../run12/restart_00010101-190000.000", 
                },                  
}
# p=11
exp_list_P11 = {
  "Dx27m_P11": { "nprcx": 10, "nprcy": 10, "Ex": 1, "Ey": 1, "Ez": 5, "porder": 11,  "analysis_elapse_time": "00:10:00", 
                  "analysis_run_no_s": 3, "analysis_run_no_e": 24, 
                  "analysis_start_time_sec": 14400, "analysis_integ_time_per_run": 7200.0,  "analysis_out_tintrv": 900.0,
                  "analysis_in_bs_filebase": "../run2/restart_00010101-040000.000", 
                },   
  "Dx13m_P11": { "nprcx": 20, "nprcy": 20, "Ex": 1, "Ey": 1, "Ez": 10, "porder": 11,  "analysis_elapse_time": "00:10:00", 
                  "analysis_run_no_s": 8, "analysis_run_no_e": 24, 
                  "analysis_start_time_sec": 7*7200, "analysis_integ_time_per_run": 7200.0,  "analysis_out_tintrv": 900.0,
                  "analysis_in_bs_filebase": "../run8/restart_00010101-160000.000", 
                }, 
  "Dx6.7m_P11": { "nprcx": 40, "nprcy": 40,  "Ex": 1, "Ey": 1, "Ez": 20, "porder": 11, "analysis_elapse_time": "00:10:00", 
                  "analysis_run_no_s": 13, "analysis_run_no_e": 36, 
                  "analysis_start_time_sec": 19*3600, "analysis_integ_time_per_run": 7200.0,  "analysis_out_tintrv": 900.0,
                  "analysis_in_bs_filebase": "../run12/restart_00010101-190000.000", 
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
