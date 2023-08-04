import os
import math
import mkconf_sub

# p=3
exp_list_p3 = {
    "Eh10Ez6P3": {"nprc": 6, "Eh": 10, "Ez":6, "porder": 3, "dt": 120, 
                   "fz": "0.00D0, 1390.57D0, 5116.68D0, 10348.47D0, 16455.20D0, 23072.18D0, 30000.00D0", 
                   "mf_alph": "2D-1", "mf_ordh": 16, "mf_alpv": "5D-2", "mf_ordv": 16, 
                   "initgp_porder": 7, "elapse_time": "00:10:00",
                   "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 }, 
    "Eh20Ez6P3": {"nprc": 24, "Eh": 20, "Ez":6, "porder": 3, "dt": 60, 
                   "fz": "0.00D0, 1390.57D0, 5116.68D0, 10348.47D0, 16455.20D0, 23072.18D0, 30000.00D0", 
                   "mf_alph": "2D-1", "mf_ordh": 16, "mf_alpv": "5D-2", "mf_ordv": 16, 
                   "initgp_porder": 7, "elapse_time": "00:20:00",
                   "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },     
    "Eh40Ez6P3": {"nprc": 96, "Eh": 40, "Ez":6, "porder": 3, "dt": 30, 
                   "fz": "0.00D0, 1390.57D0, 5116.68D0, 10348.47D0, 16455.20D0, 23072.18D0, 30000.00D0", 
                   "mf_alph": "2D-1", "mf_ordh": 16, "mf_alpv": "5D-2", "mf_ordv": 16, 
                   "initgp_porder": 7, "elapse_time": "00:40:00",
                   "regrid_nprcx": 16, "regrid_Ex": 16, "regrid_nprcy": 16, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },         
    "Eh80Ez6P3": {"nprc": 384, "Eh": 80, "Ez":6, "porder": 3, "dt": 15, 
                   "fz": "0.00D0, 1390.57D0, 5116.68D0, 10348.47D0, 16455.20D0, 23072.18D0, 30000.00D0", 
                   "mf_alph": "2D-1", "mf_ordh": 16, "mf_alpv": "5D-2", "mf_ordv": 16, 
                   "initgp_porder": 7, "elapse_time": "00:40:00",
                   "regrid_nprcx": 16, "regrid_Ex": 32, "regrid_nprcy": 16, "regrid_Ey": 16, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },             
}
# p = 7
exp_list_p7 = {
    "Eh5Ez3P7": {"nprc": 6, "Eh": 5, "Ez":3, "porder": 7, "dt": 120, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0", 
                   "mf_alph": "1D0", "mf_ordh": 16, "mf_alpv": "1D0", "mf_ordv": 16, 
                   "initgp_porder": 11, "elapse_time": "00:10:00",
                   "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 }, 
    "Eh10Ez3P7": {"nprc": 24, "Eh": 10, "Ez":3, "porder": 7, "dt": 60, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0", 
                   "mf_alph": "1D0", "mf_ordh": 16, "mf_alpv": "1D0", "mf_ordv": 16, 
                   "initgp_porder": 11, "elapse_time": "00:20:00",
                   "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },     
    "Eh20Ez3P7": {"nprc": 96, "Eh": 20, "Ez":3, "porder": 7, "dt": 30, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0", 
                   "mf_alph": "1D0", "mf_ordh": 16, "mf_alpv": "1D0", "mf_ordv": 16, 
                   "initgp_porder": 11, "elapse_time": "00:40:00",
                   "regrid_nprcx": 16, "regrid_Ex": 16, "regrid_nprcy": 16, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },         
    "Eh40Ez3P7": {"nprc": 384, "Eh": 40, "Ez":3, "porder": 7, "dt": 15, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0", 
                   "mf_alph": "1D0", "mf_ordh": 16, "mf_alpv": "1D0", "mf_ordv": 16, 
                   "initgp_porder": 11, "elapse_time": "01:30:00",
                   "regrid_nprcx": 16, "regrid_Ex": 24, "regrid_nprcy": 16, "regrid_Ey": 12, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },             
    "Eh80Ez3P7": {"nprc": 1536, "Eh": 80, "Ez":3, "porder": 7, "dt": 7.5, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0", 
                   "mf_alph": "1D0", "mf_ordh": 16, "mf_alpv": "1D0", "mf_ordv": 16, 
                   "initgp_porder": 11, "elapse_time": "03:00:00",
                   "regrid_nprcx": 16, "regrid_Ex": 32, "regrid_nprcy": 16, "regrid_Ey": 16, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },                 
}
# p = 11
exp_list_p11 = {
    "Eh4Ez3P11": {"nprc": 24, "Eh": 4, "Ez":3, "porder": 11, "dt": 60, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0", 
                   "mf_alph": "5D2", "mf_ordh": 32, "mf_alpv": "1D1", "mf_ordv": 32, 
                   "initgp_porder": 13, "elapse_time": "00:30:00",
                   "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 }, 
    "Eh8Ez3P11": {"nprc": 96, "Eh": 8, "Ez":3, "porder": 11, "dt": 30, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0", 
                   "mf_alph": "5D2", "mf_ordh": 32, "mf_alpv": "1D1", "mf_ordv": 32, 
                   "initgp_porder": 13, "elapse_time": "00:60:00",
                   "regrid_nprcx": 16, "regrid_Ex": 8, "regrid_nprcy": 16, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },     
    "Eh16Ez3P11": {"nprc": 384, "Eh": 16, "Ez":3, "porder": 11, "dt": 15, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0", 
                   "mf_alph": "5D2", "mf_ordh": 32, "mf_alpv": "1D1", "mf_ordv": 32, 
                   "initgp_porder": 13, "elapse_time": "02:00:00",
                   "regrid_nprcx": 16, "regrid_Ex": 16, "regrid_nprcy": 16, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },         
    "Eh32Ez3P11": {"nprc": 1536, "Eh": 32, "Ez":3, "porder": 11, "dt": 7.5, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0", 
                   "mf_alph": "5D2", "mf_ordh": 32, "mf_alpv": "1D1", "mf_ordv": 32, 
                   "initgp_porder": 13, "elapse_time": "04:00:00",
                   "regrid_nprcx": 16, "regrid_Ex": 16, "regrid_nprcy": 16, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
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