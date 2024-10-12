import mkconf_sub

# p=3  (b=3)
exp_list_p3 = {
    "Eh12Ez8P3": {"nprc": 24, "Eh": 12, "Ez":8, "porder": 3, "dt": 600, "dt_dyn": 75, 
                   "fz": "0.00D0, 695.07D0, 2691.74D0, 5772.72D0, 9686.27D0, 14211.85D0, 19180.79D0, 24471.90D0, 30000.00D0", 
                   "mf_alph": "1D-1", "mf_ordh": 8, "mf_alpv": "5D-2", "mf_ordv": 8, 
                   "mf_alph_ini": "5D-1", "mf_alpv_ini": "5D-2", 
                   "day_per_run": 250, "elapse_time": "06:00:00",
                   "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                   "spectra_Mt": 64, "spectra_nprc":6, "spectra_elapse_time": "00:10:00",                     
                 },
    "Eh24Ez16P3": {"nprc": 96, "Eh": 24, "Ez":16, "porder": 3, "dt": 600, "dt_dyn": 37.5, 
                   "fz": "0.00D0, 175.27D0, 695.07D0, 1542.38D0, 2691.74D0, 4112.64D0, 5772.72D0, 7640.36D0, 9686.27D0, 11884.33D0, 14211.85D0, 16649.46D0, 19180.79D0, 21792.10D0, 24471.90D0, 27210.55D0, 30000.00D0", 
                   "mf_alph": "1D-1", "mf_ordh": 8, "mf_alpv": "5D-2", "mf_ordv": 8, 
                   "mf_alph_ini": "5D-1", "mf_alpv_ini": "5D-2", 
                   "day_per_run": 125, "elapse_time": "06:00:00",
                   "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                   "spectra_Mt": 128, "spectra_nprc":24, "spectra_elapse_time": "00:20:00",  
                 },    
}
# p = 7
exp_list_p7 = {
    "Eh6Ez4P7": {"nprc": 24, "Eh": 6, "Ez":4, "porder": 7, "dt": 600, "dt_dyn": 75, 
                   "fz": "0.0000D3, 2.6917D3, 9.6863D3, 19.1808D3, 30.0000D3",
                   "mf_alph": "4D0", "mf_ordh": 16, "mf_alpv": "4D0", "mf_ordv": 16, 
                   "mf_alph_ini": "5D0", "mf_alpv_ini": "5D0", 
                   "day_per_run": 250, "elapse_time": "06:00:00",
                   "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                   "spectra_Mt": 64, "spectra_nprc":6, "spectra_elapse_time": "00:10:00",                                        
                 },  # (b=3)                  
    "Eh12Ez8P7": {"nprc": 96, "Eh": 12, "Ez":8, "porder": 7, "dt": 600, "dt_dyn": 37.5,  
                   "fz": "0.00D0, 695.07D0, 2691.74D0, 5772.72D0, 9686.27D0, 14211.85D0, 19180.79D0, 24471.90D0, 30000.00D0", 
                   "mf_alph": "4D0", "mf_ordh": 16, "mf_alpv": "4D0", "mf_ordv": 16, 
                   "mf_alph_ini": "5D0", "mf_alpv_ini": "5D0", 
                   "day_per_run": 125, "elapse_time": "06:00:00",
                   "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                   "spectra_Mt": 128, "spectra_nprc":24, "spectra_elapse_time": "00:20:00",                     
                 }, # (b=3)                     
    "Eh24Ez16P7": {"nprc": 384, "Eh": 24, "Ez":16, "porder": 7, "dt": 300, "dt_dyn": 18.75, 
                   "fz": "0.00D0, 175.27D0, 695.07D0, 1542.38D0, 2691.74D0, 4112.64D0, 5772.72D0, 7640.36D0, 9686.27D0, 11884.33D0, 14211.85D0, 16649.46D0, 19180.79D0, 21792.10D0, 24471.90D0, 27210.55D0, 30000.00D0", 
                   "mf_alph": "5D0", "mf_ordh": 16, "mf_alpv": "5D0", "mf_ordv": 16, 
                   "mf_alph_ini": "5D0", "mf_alpv_ini": "5D0", 
                   "day_per_run": 50, "elapse_time": "08:00:00",
                   "regrid_nprcx": 16, "regrid_Ex": 16, "regrid_nprcy": 16, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                   "spectra_Mt": 256, "spectra_nprc":96, "spectra_elapse_time": "00:20:00",                     
                 }, # (b=3)                 
}
# p = 11
exp_list_p11 = {
    "Eh4Ez3P11": {"nprc": 96, "Eh": 4, "Ez":3, "porder": 11, "dt": 540, "dt_dyn": 45, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0",
                   "mf_alph": "5D0", "mf_ordh": 16, "mf_alpv": "4D0", "mf_ordv": 16, 
                   "mf_alph_ini": "5D0", "mf_alpv_ini": "5D0", 
                   "day_per_run": 250, "elapse_time": "6:00:00",
                   "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                   "spectra_Mt": 64, "spectra_nprc":6, "spectra_elapse_time": "00:10:00",                                        
                 },       
    "Eh8Ez6P11": {"nprc": 384, "Eh": 8, "Ez":6, "porder": 11, "dt": 540, "dt_dyn": 22.5, 
                   "fz": "0.0D0, 1390.57D0, 5116.68D0, 10348.47D0, 16455.2D0, 23072.18D0, 30000D0",
                   "mf_alph": "4D0", "mf_ordh": 16, "mf_alpv": "4D0", "mf_ordv": 16, 
                   "mf_alph_ini": "5D0", "mf_alpv_ini": "5D0", 
                   "day_per_run": 125, "elapse_time": "6:00:00",
                   "regrid_nprcx": 16, "regrid_Ex": 8, "regrid_nprcy": 16, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:30:00",
                   "spectra_Mt": 128, "spectra_nprc":24, "spectra_elapse_time": "00:20:00",                     
                 },           
}

exp_list = {  
           **exp_list_p3, 
           **exp_list_p7, 
           **exp_list_p11,  
            }
  
#---------------------------------
for exp_name, exp_info in exp_list.items():
  mkconf_sub.mk_conf_sh( exp_name, exp_info )
