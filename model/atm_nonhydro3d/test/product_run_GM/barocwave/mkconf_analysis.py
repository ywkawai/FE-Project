import mkconf_analysis_sub

# p=3
exp_list_p3 = {
    "Eh10Ez8P3": {"nprc": 6, "Eh": 10, "Ez": 8, "porder": 3, 
                   "fz": "0.0000D3, 1.2196D3, 4.1869D3, 7.9766D3, 12.1378D3, 16.4847D3, 20.9347D3, 25.4471D3, 30.0000D3", 
                   "regrid_elapse_time": "00:30:00", 
                   "regrid_nprc": 1536, "regrid_Eh": 160,  
                 }, 
    "Eh20Ez8P3": {"nprc": 24, "Eh": 20, "Ez": 8, "porder": 3, 
                   "fz": "0.0000D3, 1.2196D3, 4.1869D3, 7.9766D3, 12.1378D3, 16.4847D3, 20.9347D3, 25.4471D3, 30.0000D3", 
                   "regrid_elapse_time": "00:30:00", 
                   "regrid_nprc": 1536, "regrid_Eh": 160,  
                 },     
    "Eh40Ez8P3": {"nprc": 96, "Eh": 40, "Ez": 8, "porder": 3, 
                   "fz": "0.0000D3, 1.2196D3, 4.1869D3, 7.9766D3, 12.1378D3, 16.4847D3, 20.9347D3, 25.4471D3, 30.0000D3", 
                   "regrid_elapse_time": "00:30:00",
                   "regrid_nprc": 1536, "regrid_Eh": 160,  
                 },         
    "Eh80Ez8P3": {"nprc": 384, "Eh": 80, "Ez": 8, "porder": 3, 
                   "fz": "0.0000D3, 1.2196D3, 4.1869D3, 7.9766D3, 12.1378D3, 16.4847D3, 20.9347D3, 25.4471D3, 30.0000D3", 
                   "regrid_elapse_time": "00:30:00",
                   "regrid_nprc": 1536, "regrid_Eh": 160,  
                 },             
}
# p = 7
exp_list_p7 = {
    "Eh5Ez4P7": {"nprc": 6, "Eh": 5, "Ez":4, "porder": 7, 
                   "fz": "0.0000D3, 4.1869D3, 12.1378D3, 20.9347D3, 30.0000D3", 
                   "regrid_elapse_time": "00:30:00",
                   "regrid_nprc": 1536, "regrid_Eh": 80,                    
                 }, 
    "Eh10Ez4P7": {"nprc": 24, "Eh": 10, "Ez":4, "porder": 7, 
                   "fz": "0.0000D3, 4.1869D3, 12.1378D3, 20.9347D3, 30.0000D3", 
                   "regrid_elapse_time": "00:30:00",
                   "regrid_nprc": 1536, "regrid_Eh": 80,                    
                 },     
    "Eh20Ez4P7": {"nprc": 96, "Eh": 20, "Ez":4, "porder": 7, 
                   "fz": "0.0000D3, 4.1869D3, 12.1378D3, 20.9347D3, 30.0000D3", 
                   "regrid_elapse_time": "00:30:00",                   
                   "regrid_nprc": 1536, "regrid_Eh": 80,                    
                 },         
    "Eh40Ez4P7": {"nprc": 384, "Eh": 40, "Ez":4, "porder": 7, 
                   "fz": "0.0000D3, 4.1869D3, 12.1378D3, 20.9347D3, 30.0000D3", 
                   "regrid_elapse_time": "00:30:00",
                   "regrid_nprc": 1536, "regrid_Eh": 80,                                       
                 },             
}
# p = 11
exp_list_p11 = {
    "Eh3Ez3P11": {"nprc": 6, "Eh": 3, "Ez":3, "porder": 11, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0", 
                   "regrid_elapse_time": "00:30:00",
                   "regrid_nprc": 1536, "regrid_Eh": 48,                                                          
                 }, 
    "Eh6Ez3P11": {"nprc": 24, "Eh": 6, "Ez":3, "porder": 11, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0", 
                   "regrid_elapse_time": "00:30:00",
                   "regrid_nprc": 1536, "regrid_Eh": 48,                                                                             
                 },     
    "Eh12Ez3P11": {"nprc": 96, "Eh": 12, "Ez":3, "porder": 11, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0", 
                   "regrid_elapse_time": "01:00:00",
                   "regrid_nprc": 1536, "regrid_Eh": 48,                       
                 },         
    "Eh24Ez3P11": {"nprc": 384, "Eh": 24, "Ez":3, "porder": 11, 
                   "fz": "0.00D0, 5116.68D0, 16455.20D0, 30000.00D0", 
                   "regrid_elapse_time": "01:00:00",
                   "regrid_nprc": 1536, "regrid_Eh": 48,                                          
                 },             
}
exp_list = {  
            **exp_list_p3,              
            **exp_list_p7,  
            **exp_list_p11,              
            }

#---------------------------------
for exp_name, exp_info in exp_list.items():
  mkconf_analysis_sub.mk_conf_sh( exp_name, exp_info )