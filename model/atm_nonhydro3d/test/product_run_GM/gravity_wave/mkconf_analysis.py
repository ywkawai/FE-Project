import mkconf_analysis_sub

# p=1
exp_list_p1 = {
    "Eh16Ez12P1": {"nprc": 6, "Eh": 16, "Ez":12, "porder": 1, 
                   "regrid_elapse_time": "00:30:00",
                 }, 
    "Eh32Ez24P1": {"nprc": 24, "Eh": 32, "Ez":24, "porder": 1, 
                   "regrid_elapse_time": "00:30:00",
                 },     
    "Eh64Ez48P1": {"nprc": 96, "Eh": 64, "Ez":48, "porder": 1, 
                   "regrid_elapse_time": "00:60:00",
                 },           
    "Eh128Ez96P1": {"nprc": 384, "Eh": 128, "Ez":96, "porder": 1, 
                   "regrid_elapse_time": "00:60:00",
                 },               
}
# p=3
exp_list_p3 = {
    "Eh8Ez6P3": {"nprc": 6, "Eh": 8, "Ez":6, "porder": 3, 
                   "regrid_elapse_time": "00:30:00",
                 }, 
    "Eh16Ez12P3": {"nprc": 24, "Eh": 16, "Ez":12, "porder": 3, 
                   "regrid_elapse_time": "00:30:00",
                 },     
    "Eh32Ez24P3": {"nprc": 96, "Eh": 32, "Ez":24, "porder": 3, 
                   "regrid_elapse_time": "00:60:00",
                 },           
    "Eh64Ez48P3": {"nprc": 384, "Eh": 64, "Ez":48, "porder": 3, 
                   "regrid_elapse_time": "00:60:00",
                 },               
}
# p = 7
exp_list_p7 = {
    "Eh4Ez3P7": {"nprc": 6, "Eh": 4, "Ez":3, "porder": 7, 
                   "regrid_elapse_time": "00:30:00",
                 }, 
    "Eh8Ez6P7": {"nprc": 24, "Eh": 8, "Ez":6, "porder": 7, 
                   "regrid_elapse_time": "00:30:00",
                 },     
    "Eh16Ez12P7": {"nprc": 96, "Eh": 16, "Ez":12, "porder": 7, 
                   "regrid_elapse_time": "00:60:00",
                 },         
    "Eh32Ez24P7": {"nprc": 384, "Eh": 32, "Ez":24, "porder": 7, 
                   "regrid_elapse_time": "00:60:00",
                 },             
}
exp_list_p7_dtx0p5 = {
    "Eh4Ez3P7_dtx0.5": {"nprc": 6, "Eh": 4, "Ez":3, "porder": 7, 
                   "regrid_elapse_time": "00:30:00",
                 }, 
    "Eh8Ez6P7_dtx0.5": {"nprc": 24, "Eh": 8, "Ez":6, "porder": 7, 
                   "regrid_elapse_time": "00:30:00",
                 },     
    "Eh16Ez12P7_dtx0.5": {"nprc": 96, "Eh": 16, "Ez":12, "porder": 7, 
                   "regrid_elapse_time": "00:60:00",
                 },         
    "Eh32Ez24P7_dtx0.5": {"nprc": 384, "Eh": 32, "Ez":24, "porder": 7, 
                   "regrid_elapse_time": "00:60:00",
                 },             
}
exp_list_p7_dtx0p25 = {
    "Eh4Ez3P7_dtx0.25": {"nprc": 6, "Eh": 4, "Ez":3, "porder": 7, 
                   "regrid_elapse_time": "00:30:00",
                 }, 
    "Eh8Ez6P7_dtx0.25": {"nprc": 24, "Eh": 8, "Ez":6, "porder": 7, 
                   "regrid_elapse_time": "00:30:00",
                 },     
    "Eh16Ez12P7_dtx0.25": {"nprc": 96, "Eh": 16, "Ez":12, "porder": 7, 
                   "regrid_elapse_time": "00:60:00",
                 },         
    "Eh32Ez24P7_dtx0.25": {"nprc": 384, "Eh": 32, "Ez":24, "porder": 7, 
                   "regrid_elapse_time": "00:60:00",
                 },             
}
# p = 11
exp_list_p11 = {
    "Eh3Ez2P11": {"nprc": 6, "Eh": 3, "Ez":2, "porder": 11, 
                   "regrid_elapse_time": "00:30:00",
                 }, 
    "Eh6Ez4P11": {"nprc": 24, "Eh": 6, "Ez":4, "porder": 11, 
                   "regrid_elapse_time": "00:30:00",
                 },     
    "Eh12Ez8P11": {"nprc": 96, "Eh": 12, "Ez":8, "porder": 11, 
                   "regrid_elapse_time": "00:30:00",
                 },         
}
exp_list_p11_dtx0p5 = {
    "Eh3Ez2P11_dtx0.5": {"nprc": 6, "Eh": 3, "Ez":2, "porder": 11, 
                   "regrid_elapse_time": "00:30:00",
                 }, 
    "Eh6Ez4P11_dtx0.5": {"nprc": 24, "Eh": 6, "Ez":4, "porder": 11, 
                   "regrid_elapse_time": "00:30:00",
                 },     
    "Eh12Ez8P11_dtx0.5": {"nprc": 96, "Eh": 12, "Ez":8, "porder": 11, 
                   "regrid_elapse_time": "00:30:00",
                 },         
}
exp_list_p11_dtx0p25 = {
    "Eh3Ez2P11_dtx0.25": {"nprc": 6, "Eh": 3, "Ez":2, "porder": 11, 
                   "regrid_elapse_time": "00:30:00",
                 }, 
    "Eh6Ez4P11_dtx0.25": {"nprc": 24, "Eh": 6, "Ez":4, "porder": 11, 
                   "regrid_elapse_time": "00:30:00",
                 },     
    "Eh12Ez8P11_dtx0.25": {"nprc": 96, "Eh": 12, "Ez":8, "porder": 11, 
                   "regrid_elapse_time": "00:30:00",
                 },         
}

exp_list = {  
            **exp_list_p1,                          
            **exp_list_p3,              
            **exp_list_p7,  
            **exp_list_p7_dtx0p5,                         
            **exp_list_p7_dtx0p25,             
            **exp_list_p11,              
            **exp_list_p11_dtx0p5, 
            **exp_list_p11_dtx0p25,            
            }

#---------------------------------
for exp_name, exp_info in exp_list.items():
  mkconf_analysis_sub.mk_conf_sh( exp_name, exp_info )