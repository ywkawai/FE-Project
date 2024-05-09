import os
import mkgraph_numerror_conv_sub as mncs

TOP_DIR="./"
OUT_DIR="./analysis_out/h25m/num_error"

p3_list = ["Eh24Ez12P3", "Eh48Ez20P3", "Eh96Ez36P3"]
p3_dof_list = {"Eh24Ez12P3": 96, "Eh48Ez20P3": 192, "Eh96Ez36P3": 384}

p7_list = ["Eh12Ez6P7", "Eh24Ez12P7", "Eh48Ez20P7"]
p7_dof_list = {"Eh12Ez6P7":96, "Eh24Ez12P7": 192,  "Eh48Ez20P7": 384}

p11_list = ["Eh8Ez5P11", "Eh16Ez8P11", "Eh32Ez14P11"]
p11_dof_list = {"Eh8Ez5P11":96, "Eh16Ez8P11": 192, "Eh32Ez14P11": 384}

#------

numerror_list = {}
dof_list = {"P3":p3_dof_list, "P7":p7_dof_list, "P11":p11_dof_list}
exp_dir="rhot_heve"
analysis_data="analysis/NUMERROR_LOG.peall"

numerror_list[f"P3"] = mncs.get_numerror_data_list(exp_dir, p3_list, analysis_data)
numerror_list[f"P7"] = mncs.get_numerror_data_list(exp_dir, p7_list, analysis_data)
numerror_list[f"P11"] = mncs.get_numerror_data_list(exp_dir, p11_list, analysis_data)

os.makedirs(f"{OUT_DIR}", exist_ok=True)

#-- time = 7200 s 
time=7200

# U --------
varname="U"

l1_error_list = {}; l2_error_list = {}; linf_error_list = {}
for key, numerror in numerror_list.items():
  l1_error_list[key], l2_error_list[key], linf_error_list[key] = mncs.resoldep_numerrordata_sel_time(varname, time, numerror_list[key], dof_list[key])

slope_dof = l1_error_list["P7"].DOF
l1_ylim = [3e-15,2e-8]
slope_L1_list = {
#  "o2": 1.5e-2*(1.0/slope_dof)**2,
  "o3": 1.5e-2*(1.0/slope_dof)**3,    
  "o4": 1e0*(1.0/slope_dof)**4,   
#  "o6": 6e2*(1.0/slope_dof)**6,    
  "o7": 8e4*(1.0/slope_dof)**7,
  "o8_2": 1.5e6*(1.0/slope_dof)**8,  
  "o12": 1e14*(1.0/slope_dof)**12}

l2_ylim = [3e-14,6e-8]
slope_L2_list = {
  "o2": 5e-4*(1.0/slope_dof)**2,
  "o3": 4e-2*(1.0/slope_dof)**3,    
  "o4": 3e0*(1.0/slope_dof)**4,   
  "o6": 7e3*(1.0/slope_dof)**6,    
  "o8": 3e8*(1.0/slope_dof)**8,
  "o8_2": 3e7*(1.0/slope_dof)**8,  
  "o12": 1e15*(1.0/slope_dof)**12}

linf_ylim = [3e-11,6e-6]
slope_Linf_list = {
  "o2": 3.5e-2*(1.0/slope_dof)**2,
  "o3": 4e-0*(1.0/slope_dof)**3,    
  "o4": 3e2*(1.0/slope_dof)**4,   
  "o7": 1e9*(1.0/slope_dof)**7,    
  "o9": .8e13*(1.0/slope_dof)**9,
#  "o10": 1.5e15*(1.0/slope_dof)**10,  
  "o12": 6e17*(1.0/slope_dof)**12}

mncs.mkgraph_num_convergence(
  [l1_error_list,l2_error_list,linf_error_list], 
  [slope_dof, slope_dof, slope_dof],  
  [slope_L1_list, slope_L2_list, slope_Linf_list], 
  [f"L1 error", f"L2 error", f"Linf error"], 
  [l1_ylim, l2_ylim, linf_ylim], f"{OUT_DIR}/numerror_{varname}_t2hr.png" )

# W --------
varname="W"

l1_error_list = {}; l2_error_list = {}; linf_error_list = {}
for key, numerror in numerror_list.items():
  l1_error_list[key], l2_error_list[key], linf_error_list[key] = mncs.resoldep_numerrordata_sel_time(varname, time, numerror_list[key], dof_list[key])

slope_dof = l1_error_list["P7"].DOF
l1_ylim = [2e-11,3e-4]
slope_L1_list = {
#  "o2": 2e0*(1.0/slope_dof)**2,
  "o3": 1.2e2*(1.0/slope_dof)**3,   
  "o4": 1.2e4*(1.0/slope_dof)**4,   
  "o8": 2e11*(1.0/slope_dof)**8,
  "o10": 2e15*(1.0/slope_dof)**10,  
  "o12": 5e18*(1.0/slope_dof)**12}

l2_ylim = [1e-9,2e-3]
slope_L2_list = {
#  "o2": 1e1*(1.0/slope_dof)**2,
  "o3": .7e3*(1.0/slope_dof)**3,     
  "o4": 7e4*(1.0/slope_dof)**4,   
#  "o4_2": 3e4*(1.0/slope_dof)**4,   
  "o8": .8e13*(1.0/slope_dof)**8,
  "o10": 5e16*(1.0/slope_dof)**10,    
  "o12": 4e19*(1.0/slope_dof)**12}

linf_ylim = [2e-6,2e-1]
slope_Linf_list = {
#  "o2": 1.5e3*(1.0/slope_dof)**2,
  "o3": 1.5e5*(1.0/slope_dof)**3,       
  "o4": 1e7*(1.0/slope_dof)**4,   
#  "o4_2": 1e6*(1.0/slope_dof)**4,   
  "o8": 1e16*(1.0/slope_dof)**8,  
  "o8_2": 1.5e15*(1.0/slope_dof)**8,
  "o12": 5e22*(1.0/slope_dof)**12}

mncs.mkgraph_num_convergence(
  [l1_error_list,l2_error_list,linf_error_list], 
  [slope_dof, slope_dof, slope_dof],  
  [slope_L1_list, slope_L2_list, slope_Linf_list], 
  [f"L1 error", f"L2 error", f"Linf error"], 
  [l1_ylim, l2_ylim, linf_ylim], f"{OUT_DIR}/numerror_{varname}_t2hr.png" )

# DDENS --------
varname="DDENS"

l1_error_list = {}; l2_error_list = {}; linf_error_list = {}
for key, numerror in numerror_list.items():
  l1_error_list[key], l2_error_list[key], linf_error_list[key] = mncs.resoldep_numerrordata_sel_time(varname, time, numerror_list[key], dof_list[key])

slope_dof = l1_error_list["P7"].DOF
l1_ylim = [1e-13,6e-7]
slope_L1_list = {
#  "o2": 4e-3*(1.0/slope_dof)**2,
  "o3": 4e-1*(1.0/slope_dof)**3,  
  "o4": 3e1*(1.0/slope_dof)**4,   
  "o7": 5e6*(1.0/slope_dof)**7,    
  "o8_2": 1e8*(1.0/slope_dof)**8,
  "o12": .7e16*(1.0/slope_dof)**12}

l2_ylim = [4e-12,3e-6]
slope_L2_list = {
#  "o2": 1.5e-2*(1.0/slope_dof)**2,
  "o3": 1.5e0*(1.0/slope_dof)**3,    
  "o4": 1.2e2*(1.0/slope_dof)**4,   
  "o5": 2e3*(1.0/slope_dof)**5,  
  "o6": 2e5*(1.0/slope_dof)**6,    
  "o8": 2e9*(1.0/slope_dof)**8,
  "o12": 5e16*(1.0/slope_dof)**12}

linf_ylim = [1e-9,4e-4]
slope_Linf_list = {
#  "o2": 5e-0*(1.0/slope_dof)**2,
  "o3": 2.5e2*(1.0/slope_dof)**3,      
  "o4": 2e4*(1.0/slope_dof)**4,   
  "o6": 2e8*(1.0/slope_dof)**6,  
#  "o7": 5e9*(1.0/slope_dof)**7,  
  "o8": 8e11*(1.0/slope_dof)**8,
  "o12": 2e19*(1.0/slope_dof)**12}

mncs.mkgraph_num_convergence(
  [l1_error_list,l2_error_list,linf_error_list], 
  [slope_dof, slope_dof, slope_dof],  
  [slope_L1_list, slope_L2_list, slope_Linf_list], 
  [f"L1 error", f"L2 error", f"Linf error"], 
  [l1_ylim, l2_ylim, linf_ylim], f"{OUT_DIR}/numerror_{varname}_t2hr.png" )

# DRHOT --------
varname="DRHOT"

l1_error_list = {}; l2_error_list = {}; linf_error_list = {}
for key, numerror in numerror_list.items():
  l1_error_list[key], l2_error_list[key], linf_error_list[key] = mncs.resoldep_numerrordata_sel_time(varname, time, numerror_list[key], dof_list[key])

slope_dof = l1_error_list["P7"].DOF
l1_ylim = [1e-12,1e-5]
slope_L1_list = {
#  "o2": 4e-3*(1.0/slope_dof)**2,
  "o3": 7e0*(1.0/slope_dof)**3,  
  "o4": 6e2*(1.0/slope_dof)**4,   
  "o7": 9e7*(1.0/slope_dof)**7,    
  "o8": 6e9*(1.0/slope_dof)**8,
  "o10": .8e14*(1.0/slope_dof)**10,    
  "o12": 1e17*(1.0/slope_dof)**12}

l2_ylim = [1e-11,1e-4]
slope_L2_list = {
#  "o2": 4e-3*(1.0/slope_dof)**2,
  "o3": 6e1*(1.0/slope_dof)**3,  
  "o4": 6e3*(1.0/slope_dof)**4,   
#  "o7": 2e9*(1.0/slope_dof)**7,    
  "o8": 7e11*(1.0/slope_dof)**8,
  "o10": .5e16*(1.0/slope_dof)**10,  
  "o12": 2e18*(1.0/slope_dof)**12}

linf_ylim = [4e-7,4e-2]
slope_Linf_list = {
#  "o2": 4e-3*(1.0/slope_dof)**2,
  "o3": 3e4*(1.0/slope_dof)**3,  
  "o4": 2e6*(1.0/slope_dof)**4,   
#  "o7": 2e9*(1.0/slope_dof)**7,    
  "o8": 1.5e15*(1.0/slope_dof)**8,
  "o9": .8e17*(1.0/slope_dof)**9,  
  "o12": 3e21*(1.0/slope_dof)**12}


mncs.mkgraph_num_convergence(
  [l1_error_list,l2_error_list,linf_error_list], 
  [slope_dof, slope_dof, slope_dof],  
  [slope_L1_list, slope_L2_list, slope_Linf_list], 
  [f"L1 error", f"L2 error", f"Linf error"], 
  [l1_ylim, l2_ylim, linf_ylim], f"{OUT_DIR}/numerror_{varname}_t2hr.png" )
