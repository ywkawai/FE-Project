import mkconf_sub

alph_deg_list = [0, 45, 90]

# p = 1
exp_list_p1 = {
    "Eh16Ez12P1": {"nprc": 6, "Eh": 16, "Ez":12, "porder": 1, "dt": 600, "initgp_porder": 3, "elapse_time": "00:20:00",
                 "regrid_nprcx": 4, "regrid_Ex": 8, "regrid_nprcy": 4, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 }, 
    "Eh32Ez12P1": {"nprc": 24, "Eh": 32, "Ez":12, "porder": 1, "dt": 300, "initgp_porder": 3, "elapse_time": "00:30:00",
                 "regrid_nprcx": 4, "regrid_Ex": 16, "regrid_nprcy": 4, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:30:00",
                 },
    "Eh64Ez12P1": {"nprc": 96, "Eh": 64, "Ez":12, "porder": 1, "dt": 150, "initgp_porder": 3, "elapse_time": "00:30:00",
                 "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:30:00",
                 },
    "Eh128Ez12P1": {"nprc": 384, "Eh": 128, "Ez":12, "porder": 1, "dt": 75, "initgp_porder": 3, "elapse_time": "00:60:00",
                 "regrid_nprcx": 8, "regrid_Ex": 32, "regrid_nprcy": 8, "regrid_Ey": 16, "regrid_porder": 3, "regrid_elapse_time": "00:30:00",
                 },            
}

# p = 3
exp_list_p3 = {
    "Eh8Ez6P3": {"nprc": 6, "Eh": 8, "Ez":6, "porder": 3, "dt": 600, "initgp_porder": 7, "elapse_time": "00:20:00",
                 "regrid_nprcx": 4, "regrid_Ex": 8, "regrid_nprcy": 4, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 }, 
    "Eh16Ez6P3": {"nprc": 24, "Eh": 16, "Ez":6, "porder": 3, "dt": 300, "initgp_porder": 7, "elapse_time": "00:30:00",
                 "regrid_nprcx": 4, "regrid_Ex": 16, "regrid_nprcy": 4, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:30:00",
                 },
    "Eh32Ez6P3": {"nprc": 96, "Eh": 32, "Ez":6, "porder": 3, "dt": 150, "initgp_porder": 7, "elapse_time": "00:30:00",
                 "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:30:00",
                 },
    "Eh64Ez6P3": {"nprc": 384, "Eh": 64, "Ez":6, "porder": 3, "dt": 75, "initgp_porder": 7, "elapse_time": "00:60:00",
                 "regrid_nprcx": 8, "regrid_Ex": 32, "regrid_nprcy": 8, "regrid_Ey": 16, "regrid_porder": 3, "regrid_elapse_time": "00:30:00",
                 },            
}

# p = 7
exp_list_p7 = {
    "Eh4Ez3P7": {"nprc": 6, "Eh": 4, "Ez":3, "porder": 7, "dt": 600, "initgp_porder": 11, "elapse_time": "00:20:00",
                 "regrid_nprcx": 4, "regrid_Ex": 8, "regrid_nprcy": 4, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 }, 
    "Eh8Ez3P7": {"nprc": 24, "Eh": 8, "Ez":3, "porder": 7, "dt": 300, "initgp_porder": 11, "elapse_time": "00:30:00",
                 "regrid_nprcx": 4, "regrid_Ex": 16, "regrid_nprcy": 4, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:30:00",
                 },
    "Eh16Ez3P7": {"nprc": 96, "Eh": 16, "Ez":3, "porder": 7, "dt": 150, "initgp_porder": 11, "elapse_time": "00:30:00",
                 "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:30:00",
                 },
    "Eh32Ez3P7": {"nprc": 384, "Eh": 32, "Ez":3, "porder": 7, "dt": 75, "initgp_porder": 11, "elapse_time": "00:60:00",
                 "regrid_nprcx": 8, "regrid_Ex": 32, "regrid_nprcy": 8, "regrid_Ey": 16, "regrid_porder": 3, "regrid_elapse_time": "00:30:00",
                 },            
}

# p = 11
exp_list_p11 = {
    "Eh2Ez3P11": {"nprc": 6, "Eh": 2, "Ez":3, "porder": 11, "dt": 600, "initgp_porder": 13, "elapse_time": "00:20:00",
                 "regrid_nprcx": 4, "regrid_Ex": 8, "regrid_nprcy": 4, "regrid_Ey": 4, "regrid_porder": 3, "regrid_elapse_time": "00:20:00",
                 }, 
    "Eh4Ez3P11": {"nprc": 24, "Eh": 4, "Ez":3, "porder": 11, "dt": 300, "initgp_porder": 13, "elapse_time": "00:30:00",
                 "regrid_nprcx": 4, "regrid_Ex": 16, "regrid_nprcy": 4, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:30:00",
                 },
    "Eh8Ez3P11": {"nprc": 96, "Eh": 8, "Ez":3, "porder": 11, "dt": 150, "initgp_porder": 13, "elapse_time": "00:30:00",
                 "regrid_nprcx": 8, "regrid_Ex": 16, "regrid_nprcy": 8, "regrid_Ey": 8, "regrid_porder": 3, "regrid_elapse_time": "00:30:00",
                 },
    "Eh16Ez3P11": {"nprc": 384, "Eh": 16, "Ez":3, "porder": 11, "dt": 75, "initgp_porder": 13, "elapse_time": "00:60:00",
                 "regrid_nprcx": 8, "regrid_Ex": 32, "regrid_nprcy": 8, "regrid_Ey": 16, "regrid_porder": 3, "regrid_elapse_time": "00:30:00",
                 },            
    "Eh16Ez3P11_dtx0.5": {"nprc": 384, "Eh": 16, "Ez":3, "porder": 11, "dt": 37.5, "initgp_porder": 13, "elapse_time": "00:60:00",
                 "regrid_nprcx": 8, "regrid_Ex": 32, "regrid_nprcy": 8, "regrid_Ey": 16, "regrid_porder": 3, "regrid_elapse_time": "00:30:00",
                 },                
}

exp_list = { **exp_list_p1, **exp_list_p3, **exp_list_p7, **exp_list_p11 }

#---------------------------------
   
for alph in alph_deg_list:
  for exp_name, exp_info in exp_list.items():
    mkconf_sub.mk_conf_jobsh( exp_name, exp_info, alph )