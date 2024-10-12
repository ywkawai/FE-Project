import mkconf_sub_hi

# p = 3 (b=3)
exp_list_p3 = {
    "Eh48Ez32P3": {"nprc": 1536, "Eh": 48, "Ez":32, "porder": 3, "dt": 300, "dt_dyn": 18.75, 
                   "fz": "0.000D3, 0.044D3, 0.175D3, 0.393D3, 0.695D3, 1.079D3, 1.542D3, 2.081D3, 2.692D3, 3.370D3, 4.113D3, 4.915D3, 5.773D3, 6.682D3, 7.640D3, 8.643D3, 9.686D3, 10.768D3, 11.884D3, 13.033D3, 14.212D3, 15.418D3, 16.649D3, 17.904D3, 19.181D3, 20.477D3, 21.792D3, 23.124D3, 24.472D3, 25.834D3, 27.211D3, 28.599D3, 30.000D3", 
                   "mf_alph": "5D0", "mf_ordh": 16, "mf_alpv": "5D0", "mf_ordv": 16, 
                   "mf_alph_ini": "5D0", "mf_alpv_ini": "5D0",  
                   "rg_in_basename": "../../Eh24Ez16P7/spinup4/restart_00010720-000000.000", 
                   "rg_out_basename": "restart_00010720-000000.000",                    
                   "rg_in_nprc": 384, "rg_in_Eh": 24, "rg_in_Ez":16, "rg_in_porder": 7,  
                   "rg_in_fz": "0.00D0, 175.27D0, 695.07D0, 1542.38D0, 2691.74D0, 4112.64D0, 5772.72D0, 7640.36D0, 9686.27D0, 11884.33D0, 14211.85D0, 16649.46D0, 19180.79D0, 21792.10D0, 24471.90D0, 27210.55D0, 30000.00D0",
                   "ini_day": 200, "integ_day": 1000, "day_per_run": 50, "elapse_time": "3:30:00", 
                   "regrid_nprcx": 16, "regrid_Ex": 16, "regrid_nprcy": 16, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                   "spectra_Mt": 512, "spectra_nprc":384, "spectra_elapse_time": "00:30:00", 
                 },                               
}

# p = 7 (b=3)
exp_list_p7 = {
    "Eh48Ez32P7": {"nprc": 1536, "Eh": 48, "Ez":32, "porder": 7, "dt": 300, "dt_dyn": 9.375, 
                   "fz": "0.000D3, 0.044D3, 0.175D3, 0.393D3, 0.695D3, 1.079D3, 1.542D3, 2.081D3, 2.692D3, 3.370D3, 4.113D3, 4.915D3, 5.773D3, 6.682D3, 7.640D3, 8.643D3, 9.686D3, 10.768D3, 11.884D3, 13.033D3, 14.212D3, 15.418D3, 16.649D3, 17.904D3, 19.181D3, 20.477D3, 21.792D3, 23.124D3, 24.472D3, 25.834D3, 27.211D3, 28.599D3, 30.000D3", 
                   "mf_alph": "5D0", "mf_ordh": 16, "mf_alpv": "5D0", "mf_ordv": 16, 
                   "mf_alph_ini": "5D0", "mf_alpv_ini": "5D0",  
                   "rg_in_basename": "../../Eh24Ez16P7/spinup4/restart_00010720-000000.000", 
                   "rg_out_basename": "restart_00010720-000000.000",                    
                   "rg_in_nprc": 384, "rg_in_Eh": 24, "rg_in_Ez":16, "rg_in_porder": 7,  
                   "rg_in_fz": "0.00D0, 175.27D0, 695.07D0, 1542.38D0, 2691.74D0, 4112.64D0, 5772.72D0, 7640.36D0, 9686.27D0, 11884.33D0, 14211.85D0, 16649.46D0, 19180.79D0, 21792.10D0, 24471.90D0, 27210.55D0, 30000.00D0",                    
                   "ini_day": 200, "integ_day": 1000, "day_per_run": 10, "elapse_time": "07:00:00", 
                   "regrid_nprcx": 16, "regrid_Ex": 16, "regrid_nprcy": 16, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                   "spectra_Mt": 512, "spectra_nprc":384, "spectra_elapse_time": "00:30:00", 
                 },                           
}

# p = 11
exp_list_p11 = {
    "Eh16Ez12P11": {"nprc": 1536, "Eh": 16, "Ez":12, "porder": 11, "dt": 540, "dt_dyn": 11.25, 
                   "fz": "0.000D3, 0.356D3, 1.391D3, 3.014D3, 5.117D3, 7.591D3, 10.348D3, 13.320D3, 16.455D3, 19.715D3, 23.072D3, 26.506D3, 30.000D3",
                   "mf_alph": "4D0", "mf_ordh": 16, "mf_alpv": "4D0", "mf_ordv": 16, 
                   "mf_alph_ini": "5D0", "mf_alpv_ini": "5D0", 
                   "rg_in_basename": "../../Eh8Ez6P11/spinup2/restart_00010720-000000.000", 
                   "rg_out_basename": "restart_00010720-000000.000",                    
                   "rg_in_nprc": 384, "rg_in_Eh": 8, "rg_in_Ez":6, "rg_in_porder": 11,  
                   "rg_in_fz": "0.0D0, 1390.57D0, 5116.68D0, 10348.47D0, 16455.2D0, 23072.18D0, 30000D0",
                   "ini_day": 200, "integ_day": 1000, "day_per_run": 50, "elapse_time": "08:00:00", 
                   "regrid_nprcx": 16, "regrid_Ex": 16, "regrid_nprcy": 16, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:40:00",
                   "spectra_Mt": 256, "spectra_nprc":96, "spectra_elapse_time": "00:30:00", 
                 },                           
}

exp_list = {  
            **exp_list_p3, 
            **exp_list_p7, 
            **exp_list_p11,  
            }
  
#---------------------------------
for exp_name, exp_info in exp_list.items():
  mkconf_sub_hi.mk_conf_sh( exp_name, exp_info )
