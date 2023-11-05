import mkconf_sub

#---
H0 = 25 # mountain height

# p = 3
exp_list_p3 = {
    "Eh24Ez12P3_h25m": {"nprc": 24, "Eh": 24, "Ez":12, "porder": 3, "dt": 0.375, 
                   "fz": "0D3, 2.000D3, 4.000D3, 6.000D3, 8.000D3, 10.000D3, 12.000D3, 14.000D3, 16.322D3, 19.018D3, 22.148D3, 25.781D3, 30.000D3", # dz_inc_ratio = 1.161
                   "mf_alph": "0.05D0", "mf_ordh": 16, "mf_alpv": "0.05D0", "mf_ordv": 16, 
                   "initgp_porder": 7, "elapse_time": "00:20:00",
                   "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 }, 
    "Eh48Ez20P3_h25m": {"nprc": 96, "Eh": 48, "Ez":20, "porder": 3, "dt": 0.1875, 
                   "fz": "0D3, 1.000D3, 2.000D3, 3.000D3, 4.000D3, 5.000D3, 6.000D3, 7.000D3, 8.000D3, 9.000D3, 10.000D3, 11.000D3, 12.000D3, 13.179D3, 14.569D3, 16.208D3, 18.139D3, 20.417D3, 23.102D3, 26.268D3, 30.000D3", # dz_inc_ratio = 1.17895
                   "mf_alph": "0.05D0", "mf_ordh": 16, "mf_alpv": "0.05D0", "mf_ordv": 16, 
                   "initgp_porder": 7, "elapse_time": "01:30:00",
                   "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },    
    "Eh96Ez36P3_h25m": {"nprc": 384, "Eh": 96, "Ez":36, "porder": 3, "dt": 0.09375, 
                   "fz": "0D3, 0.500D3, 1.000D3, 1.500D3, 2.000D3, 2.500D3, 3.000D3, 3.500D3, 4.000D3, 4.500D3, 5.000D3, 5.500D3, 6.000D3, 6.500D3, 7.000D3, 7.500D3, 8.000D3, 8.500D3, 9.000D3, 9.500D3, 10.000D3, 10.500D3, 11.000D3, 11.500D3, 12.000D3, 12.500D3, 13.092D3, 13.794D3, 14.626D3, 15.611D3, 16.779D3, 18.162D3, 19.801D3, 21.743D3, 24.044D3, 26.771D3, 30.00D3", # dz_inc_ratio = 1.18484
                   "mf_alph": "0.05D0", "mf_ordh": 16, "mf_alpv": "0.05D0", "mf_ordv": 16, 
                   "initgp_porder": 11, "elapse_time": "06:00:00",
                   "regrid_nprcx": 16, "regrid_Ex": 16, "regrid_nprcy": 16, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },             
}
# p = 7
exp_list_p7 = {
    "Eh12Ez6P7_h25m": {"nprc": 24, "Eh": 12, "Ez":6, "porder": 7, "dt": 0.375, 
                   "fz": "0D3, 4.000D3, 8.000D3, 12.000D3, 16.868D3, 22.792D3, 30.002D3",  # dz_inc_ratio = 1.217 
                   "mf_alph": "0.05D0", "mf_ordh": 16, "mf_alpv": "0.05D0", "mf_ordv": 16, 
                   "initgp_porder": 11, "elapse_time": "00:30:00",
                   "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 }, 
    "Eh24Ez12P7_h25m": {"nprc": 96, "Eh": 24, "Ez":12, "porder": 7, "dt": 0.1875, 
                   "fz": "0D3, 2.000D3, 4.000D3, 6.000D3, 8.000D3, 10.000D3, 12.000D3, 14.000D3, 16.322D3, 19.018D3, 22.148D3, 25.781D3, 30.000D3", # dz_inc_ratio = 1.161
                   "mf_alph": "0.05D0", "mf_ordh": 16, "mf_alpv": "0.05D0", "mf_ordv": 16, 
                   "initgp_porder": 11, "elapse_time": "02:00:00",
                   "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },
    "Eh48Ez20P7_h25m": {"nprc": 384, "Eh": 48, "Ez":20, "porder": 7, "dt": 0.09375, 
                   "fz": "0D3, 1.000D3, 2.000D3, 3.000D3, 4.000D3, 5.000D3, 6.000D3, 7.000D3, 8.000D3, 9.000D3, 10.000D3, 11.000D3, 12.000D3, 13.179D3, 14.569D3, 16.208D3, 18.139D3, 20.417D3, 23.102D3, 26.268D3, 30.000D3", # dz_inc_ratio = 1.17895
                   "mf_alph": "0.05D0", "mf_ordh": 16, "mf_alpv": "0.05D0", "mf_ordv": 16, 
                   "initgp_porder": 11, "elapse_time": "08:00:00",
                   "regrid_nprcx": 16, "regrid_Ex": 16, "regrid_nprcy": 16, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },         
    "Eh96Ez36P7_h25m": {"nprc": 384, "Eh": 96, "Ez":36, "porder": 7, "dt": 0.046875, 
                   "fz": "0D3, 0.500D3, 1.000D3, 1.500D3, 2.000D3, 2.500D3, 3.000D3, 3.500D3, 4.000D3, 4.500D3, 5.000D3, 5.500D3, 6.000D3, 6.500D3, 7.000D3, 7.500D3, 8.000D3, 8.500D3, 9.000D3, 9.500D3, 10.000D3, 10.500D3, 11.000D3, 11.500D3, 12.000D3, 12.500D3, 13.092D3, 13.794D3, 14.626D3, 15.611D3, 16.779D3, 18.162D3, 19.801D3, 21.743D3, 24.044D3, 26.771D3, 30.00D3", # dz_inc_ratio = 1.18484
                   "mf_alph": "0.05D0", "mf_ordh": 16, "mf_alpv": "0.05D0", "mf_ordv": 16, 
                   "initgp_porder": 11, "elapse_time": "24:00:00",
                   "regrid_nprcx": 16, "regrid_Ex": 16, "regrid_nprcy": 16, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },             
}
exp_list = {  
            **exp_list_p3,              
            **exp_list_p7,  
            }
#---------------------------------   
for exp_name, exp_info in exp_list.items():
  exp_info["h0"] = H0
  mkconf_sub.mk_conf_sh( exp_name, exp_info )