import mkconf_sub

# p = 3
exp_list_p3 = {
    "Eh24Ez10P3": {"nprc": 24, "Eh": 24, "Ez":10, "porder": 3, "dt": 0.375, 
                   "fz": "0.0D0, 2.0D3, 4.0D3, 6.0D3, 8.0D3, 10.0D3, 12.0D3, 14.0D3, 16.0D3, 18.0D3, 20.0D3",   
                   "mf_alph": "0.05D0", "mf_ordh": 16, "mf_alpv": "0.05D0", "mf_ordv": 16, 
                   "initgp_porder": 7, "elapse_time": "00:20:00",
                   "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 }, 
    "Eh48Ez20P3": {"nprc": 96, "Eh": 48, "Ez":20, "porder": 3, "dt": 0.1875, 
                   "fz": "0.0D3,  1.0D3,  2.0D3,  3.0D3,  4.0D3,  5.0D3,  6.0D3,  7.0D3,  8.0D3,  9.0D3, 10.0D3, 11.0D3, 12.0D3, 13.0D3, 14.0D3, 15.0D3, 16.0D3, 17.0D3, 18.0D3, 19.0D3, 20.0D3",   
                   "mf_alph": "0.05D0", "mf_ordh": 16, "mf_alpv": "0.05D0", "mf_ordv": 16, 
                   "initgp_porder": 7, "elapse_time": "01:30:00",
                   "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },    
}
# p = 7
exp_list_p7 = {
    "Eh12Ez5P7": {"nprc": 24, "Eh": 12, "Ez":5, "porder": 7, "dt": 0.375, 
                   "fz": "0.00D0, 4.0D3, 8.0D3, 12.0D3, 16.0D3, 20.0D3",   
                   "mf_alph": "0.05D0", "mf_ordh": 16, "mf_alpv": "0.05D0", "mf_ordv": 16, 
                   "initgp_porder": 11, "elapse_time": "00:20:00",
                   "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 }, 
    "Eh24Ez10P7": {"nprc": 96, "Eh": 24, "Ez":10, "porder": 7, "dt": 0.1875, 
                   "fz": "0.0D0, 2.0D3, 4.0D3, 6.0D3, 8.0D3, 10.0D3, 12.0D3, 14.0D3, 16.0D3, 18.0D3, 20.0D3",   
                   "mf_alph": "0.05D0", "mf_ordh": 16, "mf_alpv": "0.05D0", "mf_ordv": 16, 
                   "initgp_porder": 11, "elapse_time": "01:30:00",
                   "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },
    "Eh48Ez20P7": {"nprc": 384, "Eh": 48, "Ez":20, "porder": 7, "dt": 0.09375, 
                   "fz": "0.0D3,  1.0D3,  2.0D3,  3.0D3,  4.0D3,  5.0D3,  6.0D3,  7.0D3,  8.0D3,  9.0D3, 10.0D3, 11.0D3, 12.0D3, 13.0D3, 14.0D3, 15.0D3, 16.0D3, 17.0D3, 18.0D3, 19.0D3, 20.0D3",   
                   "mf_alph": "0.05D0", "mf_ordh": 16, "mf_alpv": "0.05D0", "mf_ordv": 16, 
                   "initgp_porder": 11, "elapse_time": "06:00:00",
                   "regrid_nprcx": 16, "regrid_Ex": 16, "regrid_nprcy": 16, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },         
}
exp_list = {  
            **exp_list_p3,              
            **exp_list_p7,  
            }
#---------------------------------   
for exp_name, exp_info in exp_list.items():
  mkconf_sub.mk_conf_sh( exp_name, exp_info )