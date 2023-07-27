import mkconf_sub

# p=3
exp_list_p3 = {
    "Eh12Ez8P3": {"nprc": 24, "Eh": 12, "Ez":8, "porder": 3, "dt": 75, 
                   "fz": "0.00D0, 695.07D0, 2691.74D0, 5772.72D0, 9686.27D0, 14211.85D0, 19180.79D0, 24471.90D0, 30000.00D0", 
                   "mf_alph": "5D-2", "mf_ordh": 8, "mf_alpv": "5D-2", "mf_ordv": 8, 
                   "mf_alph_ini": "5D-2", "mf_alpv_ini": "5D-2", 
                   "day_per_run": 250, "elapse_time": "06:00:00",
                   "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },
}
# p = 7
exp_list_p7 = {
    "Eh6Ez4P7": {"nprc": 24, "Eh": 6, "Ez":4, "porder": 7, "dt": 75, 
                   "fz": "0.0D0, 2760D0, 10340.0D0, 19528.0D0, 30000D0",
                   "mf_alph": "4D0", "mf_ordh": 16, "mf_alpv": "4D0", "mf_ordv": 16, 
                   "mf_alph_ini": "5D0", "mf_alpv_ini": "5D0", 
                   "day_per_run": 250, "elapse_time": "06:00:00",
                   "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },               
    "Eh12Ez8P7": {"nprc": 96, "Eh": 12, "Ez":8, "porder": 7, "dt": 37.5, 
                   "fz": "0.00D0, 695.07D0, 2691.74D0, 5772.72D0, 9686.27D0, 14211.85D0, 19180.79D0, 24471.90D0, 30000.00D0", 
                   "mf_alph": "4D0", "mf_ordh": 16, "mf_alpv": "4D0", "mf_ordv": 16, 
                   "mf_alph_ini": "5D0", "mf_alpv_ini": "5D0", 
                   "day_per_run": 125, "elapse_time": "06:00:00",
                   "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },                   
}
# p = 11
exp_list_p11 = {
}

exp_list = {  
            **exp_list_p3, 
            **exp_list_p7, 
            **exp_list_p11,  
            }
  
#---------------------------------
for exp_name, exp_info in exp_list.items():
  mkconf_sub.mk_conf_sh( exp_name, exp_info )
