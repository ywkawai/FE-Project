import mkconf_analysis_sub

qtrc_name_list = [
#  "cosine_bells", 
  "gaussian_hills", 
  "slotted_cylinders"
]
# p=3
exp_list_p1 = {     
    "Eh16Ez12P1": {"nprc": 6, "Eh": 16, "Ez":12, "porder": 1, "regrid_elapse_time": "00:10:00", }, 
    "Eh32Ez12P1": {"nprc": 24, "Eh": 32, "Ez":12, "porder": 1, "regrid_elapse_time": "00:10:00", },     
    "Eh64Ez12P1": {"nprc": 96, "Eh": 64, "Ez":12, "porder": 1, "regrid_elapse_time": "00:10:00", },         
    "Eh128Ez12P1": {"nprc": 384, "Eh": 128, "Ez":12, "porder": 1, "regrid_elapse_time": "00:10:00", },             
}
# p=3
exp_list_p3 = {     
    "Eh8Ez6P3": {"nprc": 6, "Eh": 8, "Ez":6, "porder": 3, "regrid_elapse_time": "00:10:00", }, 
    "Eh16Ez6P3": {"nprc": 24, "Eh": 16, "Ez":6, "porder": 3, "regrid_elapse_time": "00:10:00", },     
    "Eh32Ez6P3": {"nprc": 96, "Eh": 32, "Ez":6, "porder": 3, "regrid_elapse_time": "00:10:00", },         
    "Eh64Ez6P3": {"nprc": 384, "Eh": 64, "Ez":6, "porder": 3, "regrid_elapse_time": "00:10:00", },             
}
# p = 7
exp_list_p7 = {
    "Eh4Ez3P7": {"nprc": 6, "Eh": 4, "Ez":3, "porder": 7, "regrid_elapse_time": "00:10:00", }, 
    "Eh8Ez3P7": {"nprc": 24, "Eh": 8, "Ez":3, "porder": 7, "regrid_elapse_time": "00:10:00", },     
    "Eh16Ez3P7": {"nprc": 96, "Eh": 16, "Ez":3, "porder": 7, "regrid_elapse_time": "00:10:00", },         
    "Eh32Ez3P7": {"nprc": 384, "Eh": 32, "Ez":3, "porder": 7, "regrid_elapse_time": "00:10:00", },             
    # "Eh4Ez3P7_check": {"nprc": 6, "Eh": 4, "Ez":3, "porder": 7, "regrid_elapse_time": "00:10:00", }, 
    # "Eh8Ez3P7_check": {"nprc": 24, "Eh": 8, "Ez":3, "porder": 7, "regrid_elapse_time": "00:10:00", },     
    # "Eh16Ez3P7_check": {"nprc": 96, "Eh": 16, "Ez":3, "porder": 7, "regrid_elapse_time": "00:10:00", },         
    "Eh32Ez3P7_check": {"nprc": 384, "Eh": 32, "Ez":3, "porder": 7, "regrid_elapse_time": "00:10:00", },                 
}
exp_list_p7_MF = {
    "Eh4Ez3P7_MF": {"nprc": 6, "Eh": 4, "Ez":3, "porder": 7, "regrid_elapse_time": "00:10:00", }, 
    "Eh8Ez3P7_MF": {"nprc": 24, "Eh": 8, "Ez":3, "porder": 7, "regrid_elapse_time": "00:10:00", },     
    "Eh16Ez3P7_MF": {"nprc": 96, "Eh": 16, "Ez":3, "porder": 7, "regrid_elapse_time": "00:10:00", },         
    "Eh32Ez3P7_MF": {"nprc": 384, "Eh": 32, "Ez":3, "porder": 7, "regrid_elapse_time": "00:10:00", },             
}
# p = 11
exp_list_p11 = {
    "Eh2Ez3P11": {"nprc": 6, "Eh": 2, "Ez":3, "porder": 11, "regrid_elapse_time": "00:10:00", }, 
    "Eh4Ez3P11": {"nprc": 24, "Eh": 4, "Ez":3, "porder": 11, "regrid_elapse_time": "00:10:00", },     
    "Eh8Ez3P11": {"nprc": 96, "Eh": 8, "Ez":3, "porder": 11, "regrid_elapse_time": "00:10:00", },         
    "Eh16Ez3P11": {"nprc": 384, "Eh": 16, "Ez":3, "porder": 11, "regrid_elapse_time": "00:10:00", },             
}
exp_list_p11_MF = {
    "Eh2Ez3P11_MF": {"nprc": 6, "Eh": 2, "Ez":3, "porder": 11, "regrid_elapse_time": "00:10:00", }, 
    "Eh4Ez3P11_MF": {"nprc": 24, "Eh": 4, "Ez":3, "porder": 11, "regrid_elapse_time": "00:10:00", },     
    "Eh8Ez3P11_MF": {"nprc": 96, "Eh": 8, "Ez":3, "porder": 11, "regrid_elapse_time": "00:10:00", },         
    "Eh16Ez3P11_MF": {"nprc": 384, "Eh": 16, "Ez":3, "porder": 11, "regrid_elapse_time": "00:10:00", },             
}
exp_list = {  
            **exp_list_p1,                                                           
            **exp_list_p3,                                               
            **exp_list_p7,  
            **exp_list_p7_MF,              
            **exp_list_p11, 
            **exp_list_p11_MF,                           
            }

#---------------------------------
for qtrc_name in qtrc_name_list:
  for exp_name, exp_info in exp_list.items():
    mkconf_analysis_sub.mk_conf_sh( qtrc_name, exp_name, exp_info )