import mkconf_analysis_sub

# p=3
exp_list_p3 = {
}
# p = 7
exp_list_p7 = {
    "Eh4Ez3P7": {"nprc": 6, "Eh": 4, "Ez":3, "porder": 7, 
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