import mkconf_analysis_sub

# p=3
exp_list_p3 = {
    "Eh10Ez6P3": {"nprc": 6, "Eh": 10, "Ez": 6, "porder": 3, 
                   "fz": "0.00D0, 1390.57D0, 5116.68D0, 10348.47D0, 16455.20D0, 23072.18D0, 30000.00D0", 
                   "regrid_elapse_time": "00:30:00",
                 }, 
    "Eh20Ez6P3": {"nprc": 24, "Eh": 20, "Ez": 6, "porder": 3, 
                   "fz": "0.00D0, 1390.57D0, 5116.68D0, 10348.47D0, 16455.20D0, 23072.18D0, 30000.00D0", 
                   "regrid_elapse_time": "00:30:00",
                 },     
    "Eh40Ez6P3": {"nprc": 96, "Eh": 40, "Ez": 6, "porder": 3, 
                   "fz": "0.00D0, 1390.57D0, 5116.68D0, 10348.47D0, 16455.20D0, 23072.18D0, 30000.00D0", 
                   "regrid_elapse_time": "00:30:00",
                 },         
}
# p = 7
exp_list_p7 = {
    "Eh5Ez3P7": {"nprc": 6, "Eh": 5, "Ez":3, "porder": 7, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0", 
                   "regrid_elapse_time": "00:30:00",
                 }, 
    "Eh10Ez3P7": {"nprc": 24, "Eh": 10, "Ez":3, "porder": 7, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0", 
                   "regrid_elapse_time": "00:30:00",
                 },     
    "Eh20Ez3P7": {"nprc": 96, "Eh": 20, "Ez":3, "porder": 7, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0", 
                   "regrid_elapse_time": "00:30:00",
                 },         
}
exp_list = {  
            **exp_list_p3,              
            **exp_list_p7,  
            }

#---------------------------------
for exp_name, exp_info in exp_list.items():
  mkconf_analysis_sub.mk_conf_sh( exp_name, exp_info )