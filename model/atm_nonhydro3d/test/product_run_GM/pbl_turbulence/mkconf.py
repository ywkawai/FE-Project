import mkconf_sub

# p = 7
exp_list_p7 = {
    "Eh64Ez34P7": {"nprc": 1536, "Eh": 64, "Ez":34, "porder": 7, "dt": 0.0125, "dt_dyn": 0.0125, 
                  "fz": "0D0, 80.0D0, 160.0D0, 240.0D0, 320.0D0, 400.0D0, 480.0D0, 560.0D0, 640.0D0, 720.0D0, 800.0D0, 880.0D0, 960.0D0, 1040.0D0, 1120.0D0, 1200.0D0, 1280.0D0, 1360.0D0, 1440.0D0, 1520.0D0, 1600.0D0, 1680.0D0, 1760.0D0, 1840.0D0, 1920.0D0, 2004.3136D0, 2093.17338931D0, 2186.82449846D0, 2285.52527542D0, 2389.54799827D0, 2499.17962634D0, 2614.72259179D0, 2736.49563394D0, 2864.83467852D0, 3000.0D0", 
                  "mf_alph": 0.001, "mf_ordh": 32, "mf_alpv": 0.001, "mf_ordv": 32, 
                  "mf_alph_ini": 0.001, "mf_alpv_ini": 0.001,                   
                  "rplanet": 3400, "shallow_atm_approx": True, "elapse_time": "8:00:00",
                   "ini_hour": 0.0, "hour_per_run": 0.5, "run_num": 8, "hist_int_sec": 300, "runno_analysis": 8, "hist_int_sec_analysis": 60, 
                  "regrid_nprcx": 16, "regrid_Ex": 24, "regrid_nprcy": 8, "regrid_Ey": 24, "regrid_porder": 3, "regrid_elapse_time": "01:00:00",
                 },     
    "Eh64Ez34P7_deepatm": {"nprc": 1536, "Eh": 64, "Ez":34, "porder": 7, "dt": 0.0125, "dt_dyn": 0.0125, 
                  "fz": "0D0, 80.0D0, 160.0D0, 240.0D0, 320.0D0, 400.0D0, 480.0D0, 560.0D0, 640.0D0, 720.0D0, 800.0D0, 880.0D0, 960.0D0, 1040.0D0, 1120.0D0, 1200.0D0, 1280.0D0, 1360.0D0, 1440.0D0, 1520.0D0, 1600.0D0, 1680.0D0, 1760.0D0, 1840.0D0, 1920.0D0, 2004.3136D0, 2093.17338931D0, 2186.82449846D0, 2285.52527542D0, 2389.54799827D0, 2499.17962634D0, 2614.72259179D0, 2736.49563394D0, 2864.83467852D0, 3000.0D0", 
                  "mf_alph": 0.001, "mf_ordh": 32, "mf_alpv": 0.001, "mf_ordv": 32, 
                  "mf_alph_ini": 0.001, "mf_alpv_ini": 0.001,  
                  "rplanet": 3400, "shallow_atm_approx": False, "elapse_time": "8:00:00",
                  "ini_hour": 0.0, "hour_per_run": 0.5, "run_num": 8, "hist_int_sec": 300, "runno_analysis": 8, "hist_int_sec_analysis": 60, 
                  "regrid_nprcx": 16, "regrid_Ex": 24, "regrid_nprcy": 8, "regrid_Ey": 24, "regrid_porder": 3, "regrid_elapse_time": "01:00:00",
                 },         
}
exp_list = {  
            # **exp_list_p1,                      
            # **exp_list_p3,          
            **exp_list_p7,  
            # **exp_list_p11,              
            }
#---------------------------------   
for exp_name, exp_info in exp_list.items():
  mkconf_sub.mk_conf_sh( exp_name, exp_info )