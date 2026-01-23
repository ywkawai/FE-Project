import mkconf_sub_DNS as mkconf_sub

exp_info_common = {"dns_nu": "4D-1", "dns_mu": "5.714D-1", "ini_PT": "300.5D0",  "node_shape": "", 
                   "const_hflx": "200D0", "comp_comm_overlap": True, "StabCoef_bnd": "7.5D0" }

# p=3
exp_list_P3 = {
  "Dx6.3m_P3": { "node_shape": "", "nprcx": 32, "nprcy": 32, "Ex": 4, "Ey": 4, "Ez": 64, "porder": 3, "dt_sec": "1.25D-1", "dt_dyn_sec": "1.25D-2", 
                  "runno_s": 13,  "start_time_sec" : 19*3600,                    
                  "integ_hour": 48, "hr_per_run": 1, "elapse_time": "02:30:00",
                  "modal_filter_flag": True, "mf_ordh": 24, "mf_alph": "2D-3", "mf_ordv": 24, "mf_alpv": "2D-3", 
                  "hist_int_sec": "900D0", "initgp_porder": 11, 
                  "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 16, "regrid_porder": 3, "regrid_elapse_time": "00:15:00",                                     
                },      
}
# p=7
exp_list_P7 = {
  "Dx6.3m_P7": { "node_shape": "", "nprcx": 32, "nprcy": 32, "Ex": 2, "Ey": 2, "Ez": 32, "porder": 7, "dt_sec": "7.8125D-2", "dt_dyn_sec": "7.8125D-3", 
                  "runno_s": 8,  "start_time_sec" : 7*7200,                    
                  "integ_hour": 48, "hr_per_run": 1, "elapse_time": "02:30:00",
                  "modal_filter_flag": True, "mf_ordh": 32, "mf_alph": "7.5D-1", "mf_ordv": 32, "mf_alpv": "7.5D-1", 
                  "hist_int_sec": "900D0", "initgp_porder": 11, 
                  "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 7, "regrid_elapse_time": "00:15:00",                                     
                },      
  "Dx3.1m_P7": { "node_shape": "", "nprcx": 64, "nprcy": 64, "Ex": 2, "Ey": 2, "Ez": 64, "porder": 7, "dt_sec": "3.90625D-2", "dt_dyn_sec": "3.90625D-3", 
                  "runno_s": 13,  "start_time_sec" : 19*3600, 
                  "integ_hour": 48, "hr_per_run": 0.5, "elapse_time": "3:30:00",
                  "modal_filter_flag": True, "mf_ordh": 32, "mf_alph": "7.5D-1", "mf_ordv": 32, "mf_alpv": "7.5D-1", 
                  "hist_int_sec": "900D0", "initgp_porder": 11, 
                  "regrid_nprcx": 16, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 16, "regrid_porder": 7, "regrid_elapse_time": "00:15:00",                                     
                },      
}
# p=11
exp_list_P11 = {
  "Dx6.7m_P11": { "nprcx": 40, "nprcy": 40, "Ex": 1, "Ey": 1, "Ez": 20, "porder": 11, "dt_sec": "6.25D-2", "dt_dyn_sec": "6.25D-3", 
                  "runno_s": 13,  "start_time_sec" : 19*3600,                                                              
                  "integ_hour": 72, "hr_per_run": 2, "elapse_time": "03:30:00",
                  "modal_filter_flag": True, "mf_ordh": 40, "mf_alph": "2D1", "mf_ordv": 40, "mf_alpv": "2D1", 
                  "hist_int_sec": "900D0", "initgp_porder": 13, 
                  "regrid_nprcx": 20, "regrid_Ex": 2, "regrid_nprcy": 20, "regrid_Ey": 2, "regrid_porder": 11, "regrid_elapse_time": "00:10:00",                                                                         
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
  mkconf_sub.mk_conf_sh( exp_name, exp_info, "" )
