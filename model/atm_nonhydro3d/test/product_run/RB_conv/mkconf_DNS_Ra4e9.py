import mkconf_sub_DNS as mkconf_sub

exp_info_common = {"dns_nu": "1.512D-1", "dns_mu": "1.8D-1"}
# p=3
exp_list_P3 = {
  "Dx25m_P3": { "nprcx": 8, "nprcy": 8, "Eh": 4, "Ez": 16, "porder": 3, "dt_sec": "3D-1", "dt_dyn_sec": "3D-2", 
                  "integ_hour": 72, "hr_per_run": 2, "elapse_time": "01:00:00",
                  "modal_filter_flag": True, "mf_ordh": 32, "mf_alph": "1D-1", "mf_ordv": 32, "mf_alpv": "1D-1",
                  "hist_int_sec": 600, "initgp_porder": 7, 
                },
  "Dx12.5m_P3": { "nprcx": 16, "nprcy": 16, "Eh": 4, "Ez": 32, "porder": 3, "dt_sec": "1.5D-1", "dt_dyn_sec": "1.5D-2", 
                  "integ_hour": 72, "hr_per_run": 2, "elapse_time": "04:00:00",
                  "modal_filter_flag": True, "mf_ordh": 32, "mf_alph": "1D-1", "mf_ordv": 32, "mf_alpv": "1D-1",
                  "hist_int_sec": 600, "initgp_porder": 7, 
                },          
}
# p=7
exp_list_P7 = {
  "Dx25m_P7": { "nprcx": 8, "nprcy": 8, "Eh": 2, "Ez": 8, "porder": 7, "dt_sec": "3D-1", "dt_dyn_sec": "3D-2", 
                  "integ_hour": 72, "hr_per_run": 2, "elapse_time": "01:00:00",
                  "modal_filter_flag": True, "mf_ordh": 32, "mf_alph": "2D0", "mf_ordv": 32, "mf_alpv": "2D0",
                  "hist_int_sec": 600, "initgp_porder": 11, 
                },
  "Dx12.5m_P7": { "nprcx": 16, "nprcy": 16, "Eh": 2, "Ez": 16, "porder": 7, "dt_sec": "1.5D-1", "dt_dyn_sec": "1.5D-2", 
                  "integ_hour": 72, "hr_per_run": 2, "elapse_time": "04:00:00",
                  "modal_filter_flag": True, "mf_ordh": 32, "mf_alph": "1D0", "mf_ordv": 32, "mf_alpv": "1D0",
                  "hist_int_sec": 600, "initgp_porder": 11, 
                },          
}
exp_list = {  
            **exp_list_P3,              
            **exp_list_P7,  
            }
#---------------------------------

for exp_name, exp_info in exp_list.items():
  exp_info.update(exp_info_common)
  mkconf_sub.mk_conf_sh( exp_name, exp_info, "_Ra4e9" )
