import mkconf_sub_LES as mkconf_sub

exp_info_common = {"runno_s": 1}
# p=3
exp_list_P3 = {
  "Dx25m_P3": { "nprcx": 8, "nprcy": 8, "Ex": 4, "Ey": 4, "Ez": 16, "porder": 3, "dt_sec": "3D-1", "dt_dyn_sec": "3D-2", 
                  "integ_hour": 72, "hr_per_run": 2, "elapse_time": "01:00:00",
                  "modal_filter_flag": True, "mf_ordh": 32, "mf_alph": "1D-3", "mf_ordv": 32, "mf_alpv": "1D-3",
                  "hist_int_sec": "600D0", "initgp_porder": 7, 
                  "regrid_nprcx": 8, "regrid_Ex": 4, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:15:00", 
                },    
  "Dx12.5m_P3": { "nprcx": 16, "nprcy": 16,  "Ex": 4, "Ey": 4, "Ez": 32, "porder": 3, "dt_sec": "1.5D-1", "dt_dyn_sec": "1.5D-2", 
                  "integ_hour": 72, "hr_per_run": 2, "elapse_time": "04:00:00",
                  "modal_filter_flag": True, "mf_ordh": 32, "mf_alph": "1D-3", "mf_ordv": 32, "mf_alpv": "1D-3",
                  "hist_int_sec": "600D0", "initgp_porder": 7, 
                  "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:15:00",                   
                },          
}
# p=7
exp_list_P7 = {
  "Dx25m_P7": { "nprcx": 8, "nprcy": 8, "Ex": 2, "Ey": 2, "Ez": 8, "porder": 7, "dt_sec": "3D-1", "dt_dyn_sec": "3D-2", 
                  "integ_hour": 72, "hr_per_run": 2, "elapse_time": "01:00:00",
                  "modal_filter_flag": True, "mf_ordh": 32, "mf_alph": "1D-3", "mf_ordv": 32, "mf_alpv": "1D-3",
                  "hist_int_sec": "600D0", "initgp_porder": 11, 
                  "regrid_nprcx": 8, "regrid_Ex": 4, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:15:00",                   
                },    
  "Dx12.5m_P7": { "nprcx": 16, "nprcy": 16, "Ex": 2, "Ey": 2, "Ez": 16, "porder": 7, "dt_sec": "1.5D-1", "dt_dyn_sec": "1.5D-2", 
                  "integ_hour": 72, "hr_per_run": 2, "elapse_time": "04:00:00",
                  "modal_filter_flag": True, "mf_ordh": 32, "mf_alph": "1D-3", "mf_ordv": 32, "mf_alpv": "1D-3",
                  "hist_int_sec": "600D0", "initgp_porder": 11, 
                  "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:15:00",                                     
                },          
}
# p=11
exp_list_P11 = {
  "Dx27m_P11": { "nprcx": 10, "nprcy": 10, "Ex": 1, "Ey": 1, "Ez": 5, "porder": 11, "dt_sec": "2.5D-1", "dt_dyn_sec": "2.5D-2", 
                  "integ_hour": 72, "hr_per_run": 2, "elapse_time": "02:00:00",
                  "modal_filter_flag": True, "mf_ordh": 32, "mf_alph": "8D-4", "mf_ordv": 32, "mf_alpv": "8D-4",
                  "hist_int_sec": "600D0", "initgp_porder": 13, 
                  "regrid_nprcx": 10, "regrid_Ex": 1, "regrid_nprcy": 10, "regrid_Ey": 1, "regrid_porder": 11, "regrid_elapse_time": "00:15:00",                                     
                },    
  "Dx13m_P11": { "nprcx": 20, "nprcy": 20, "Ex": 1, "Ey": 1, "Ez": 10, "porder": 11, "dt_sec": "1.25D-1", "dt_dyn_sec": "1.25D-2", 
                  "integ_hour": 72, "hr_per_run": 2, "elapse_time": "05:00:00",
                  "modal_filter_flag": True, "mf_ordh": 32, "mf_alph": "8D-4", "mf_ordv": 32, "mf_alpv": "8D-4",
                  "hist_int_sec": "600D0", "initgp_porder": 13, 
                  "regrid_nprcx": 10, "regrid_Ex": 2, "regrid_nprcy": 10, "regrid_Ey": 2, "regrid_porder": 11, "regrid_elapse_time": "00:15:00",                                                       
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
  mkconf_sub.mk_conf_sh( exp_name, exp_info )
