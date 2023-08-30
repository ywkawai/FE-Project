import mkconf_sub

# p = 1
exp_list_p1 = {
    "Eh16Ez12P1": {"nprc": 6, "Eh": 16, "Ez":12, "porder": 1, "dt": 180, "dt_dyn": 180, 
                  "modal_filter_flag": False, "mf_alph": "-", "mf_ordh": 0, "mf_alpv": "-", "mf_ordv": 0, 
                  "initgp_porder": 3, "elapse_time": "00:20:00",
                  "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 }, 
    "Eh32Ez24P1": {"nprc": 24, "Eh": 32, "Ez":24, "porder": 1, "dt": 180, "dt_dyn": 90, 
                  "modal_filter_flag": False, "mf_alph": "-", "mf_ordh": 0, "mf_alpv": "-", "mf_ordv": 0, 
                  "initgp_porder": 3, "elapse_time": "00:20:00",
                  "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },       
    "Eh64Ez48P1": {"nprc": 96, "Eh": 64, "Ez":48, "porder": 1, "dt": 180, "dt_dyn": 45, 
                  "modal_filter_flag": False, "mf_alph": "-", "mf_ordh": 0, "mf_alpv": "-", "mf_ordv": 0, 
                  "initgp_porder": 3, "elapse_time": "01:30:00",
                  "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },           
    "Eh128Ez96P1": {"nprc": 384, "Eh": 128, "Ez":96, "porder": 1, "dt": 180, "dt_dyn": 22.5, 
                  "modal_filter_flag": False, "mf_alph": "-", "mf_ordh": 0, "mf_alpv": "-", "mf_ordv": 0, 
                  "initgp_porder": 3, "elapse_time": "06:00:00",
                  "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },               
}
# p = 3
exp_list_p3 = {
    "Eh8Ez6P3": {"nprc": 6, "Eh": 8, "Ez":6, "porder": 3, "dt": 120, "dt_dyn": 120, 
                  "modal_filter_flag": False, "mf_alph": "-", "mf_ordh": 0, "mf_alpv": "-", "mf_ordv": 0, 
                  "initgp_porder": 7, "elapse_time": "00:20:00",
                  "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 }, 
    "Eh16Ez12P3": {"nprc": 24, "Eh": 16, "Ez":12, "porder": 3, "dt": 120, "dt_dyn": 60, 
                  "modal_filter_flag": False, "mf_alph": "-", "mf_ordh": 0, "mf_alpv": "-", "mf_ordv": 0, 
                  "initgp_porder": 7, "elapse_time": "00:20:00",
                  "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },       
    "Eh32Ez24P3": {"nprc": 96, "Eh": 32, "Ez":24, "porder": 3, "dt": 120, "dt_dyn": 30, 
                  "modal_filter_flag": False, "mf_alph": "-", "mf_ordh": 0, "mf_alpv": "-", "mf_ordv": 0, 
                  "initgp_porder": 7, "elapse_time": "01:30:00",
                  "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },           
    "Eh64Ez48P3": {"nprc": 384, "Eh": 64, "Ez":48, "porder": 3, "dt": 120, "dt_dyn": 15, 
                  "modal_filter_flag": False, "mf_alph": "-", "mf_ordh": 0, "mf_alpv": "-", "mf_ordv": 0, 
                  "initgp_porder": 7, "elapse_time": "06:00:00",
                  "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },               
}
# p = 7
exp_list_p7 = {
    "Eh4Ez3P7": {"nprc": 6, "Eh": 4, "Ez":3, "porder": 7, "dt": 120, "dt_dyn": 120, 
                  "modal_filter_flag": False, "mf_alph": "-", "mf_ordh": 0, "mf_alpv": "-", "mf_ordv": 0, 
                  "initgp_porder": 11, "elapse_time": "00:20:00",
                  "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 }, 
    "Eh8Ez6P7": {"nprc": 24, "Eh": 8, "Ez":6, "porder": 7, "dt": 120, "dt_dyn": 60, 
                  "modal_filter_flag": False, "mf_alph": "-", "mf_ordh": 0, "mf_alpv": "-", "mf_ordv": 0, 
                  "initgp_porder": 11, "elapse_time": "00:20:00",
                  "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },     
    "Eh16Ez12P7": {"nprc": 96, "Eh": 16, "Ez":12, "porder": 7, "dt": 120, "dt_dyn": 30, 
                  "modal_filter_flag": False, "mf_alph": "-", "mf_ordh": 0, "mf_alpv": "-", "mf_ordv": 0, 
                  "initgp_porder": 11, "elapse_time": "01:30:00",
                  "regrid_nprcx": 16, "regrid_Ex": 8, "regrid_nprcy": 16, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },         
    "Eh32Ez24P7": {"nprc": 384, "Eh": 32, "Ez":24, "porder": 7, "dt": 120, "dt_dyn": 15, 
                  "modal_filter_flag": False, "mf_alph": "-", "mf_ordh": 0, "mf_alpv": "-", "mf_ordv": 0, 
                  "initgp_porder": 11, "elapse_time": "02:00:00",
                  "regrid_nprcx": 16, "regrid_Ex": 16, "regrid_nprcy": 16, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 }, 
    "Eh64Ez48P7": {"nprc": 1536, "Eh": 64, "Ez":48, "porder": 7, "dt": 120, "dt_dyn": 7.5, 
                  "modal_filter_flag": False, "mf_alph": "-", "mf_ordh": 0, "mf_alpv": "-", "mf_ordv": 0, 
                  "initgp_porder": 11, "elapse_time": "05:00:00",
                  "regrid_nprcx": 16, "regrid_Ex": 16, "regrid_nprcy": 16, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },                 
}
# p = 11
exp_list_p11 = {
    "Eh3Ez2P11": {"nprc": 6, "Eh": 3, "Ez":2, "porder": 11, "dt": 100, "dt_dyn": 100, 
                  "modal_filter_flag": False, "mf_alph": "-", "mf_ordh": 0, "mf_alpv": "-", "mf_ordv": 0, 
                  "initgp_porder": 13, "elapse_time": "00:20:00",
                  "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 }, 
    "Eh6Ez4P11": {"nprc": 24, "Eh": 6, "Ez":4, "porder": 11, "dt": 50, "dt_dyn": 50, 
                  "modal_filter_flag": False, "mf_alph": "-", "mf_ordh": 0, "mf_alpv": "-", "mf_ordv": 0, 
                  "initgp_porder": 13, "elapse_time": "00:40:00",
                  "regrid_nprcx": 8, "regrid_Ex": 8, "regrid_nprcy": 8, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },     
    "Eh12Ez8P11": {"nprc": 96, "Eh": 12, "Ez":8, "porder": 11, "dt": 25, "dt_dyn": 25, 
                  "modal_filter_flag": False, "mf_alph": "-", "mf_ordh": 0, "mf_alpv": "-", "mf_ordv": 0, 
                  "initgp_porder": 13, "elapse_time": "03:00:00",
                  "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 },         
}
exp_list = {  
            **exp_list_p1,                      
            **exp_list_p3,          
            **exp_list_p7,  
            **exp_list_p11,              
            }
#---------------------------------   
for exp_name, exp_info in exp_list.items():
  mkconf_sub.mk_conf_sh( exp_name, exp_info )