import os
import mkgraph_numerror_conv_sub as mncs

TOP_DIR="./"
OUT_DIR="./analysis_out/h25m/num_error"

p3_list = ["Eh96Ez15P3", "Eh192Ez30P3", "Eh384Ez60P3", "Eh768Ez120P3", "Eh1536Ez240P3"]
p3_dof_list = {"Eh96Ez15P3":384, "Eh192Ez30P3":768, "Eh384Ez60P3":1536, "Eh768Ez120P3":3072, "Eh1536Ez240P3":6144}

p3_MFoff_list = ["Eh96Ez15P3_MFoff", "Eh192Ez30P3_MFoff", "Eh384Ez60P3_MFoff", "Eh768Ez120P3_MFoff", "Eh1536Ez240P3_MFoff"]
p3_MFoff_dof_list = {"Eh96Ez15P3_MFoff":384, "Eh192Ez30P3_MFoff":768, "Eh384Ez60P3_MFoff":1536, "Eh768Ez120P3_MFoff":3072, "Eh1536Ez240P3_MFoff":6144}

p7_list = ["Eh48Ez6P7", "Eh96Ez12P7", "Eh192Ez20P7"]
p7_dof_list = {"Eh48Ez6P7":384, "Eh96Ez12P7":768, "Eh192Ez20P7":1536}

p11_list = ["Eh32Ez5P11", "Eh64Ez10P11", "Eh128Ez20P11"]
p11_dof_list = {"Eh32Ez5P11":384, "Eh64Ez10P11": 768, "Eh128Ez20P11": 1536}

#------

numerror_list = {}
dof_list = {
    "P3":p3_dof_list, 
    "P3_MFoff":p3_dof_list,     
    "P7":p7_dof_list, 
    "P11":p11_dof_list
}
exp_dir="rhot_heve"
analysis_data="analysis/NUMERROR_LOG.peall"
var_list = ["DDENS", "U", "W", "THERM", "GsqrtDDENS", "GsqrtMOMX", "GsqrtMOMZ"]

numerror_list[f"P3"] = mncs.get_numerror_data_list(exp_dir, p3_list, analysis_data, var_list)
numerror_list[f"P3_MFoff"] = mncs.get_numerror_data_list(exp_dir, p3_MFoff_list, analysis_data, var_list)
numerror_list[f"P7"] = mncs.get_numerror_data_list(exp_dir, p7_list, analysis_data, var_list )
numerror_list[f"P11"] = mncs.get_numerror_data_list(exp_dir, p11_list, analysis_data, var_list )

os.makedirs(f"{OUT_DIR}", exist_ok=True)

#-- time = 7200 s 
time=7200

# U --------
varname="U"

l1_error_list = {}; l2_error_list = {}; linf_error_list = {}
for key, numerror in numerror_list.items():
  l1_error_list[key], l2_error_list[key], linf_error_list[key] = mncs.resoldep_numerrordata_sel_time(varname, time, numerror_list[key], dof_list[key])

slope_dof = l1_error_list["P3"].DOF
l1_ylim = [1e-11,5e-4]
slope_L1_list = {
  "o2": 1e2*(1.0/slope_dof)**2,
  "o3_2": 1.5e4*(1.0/slope_dof)**3,     
  "o4": 4e6*(1.0/slope_dof)**4,   
  "o8": 4e15*(1.0/slope_dof)**8,
  "o12": 3e25*(1.0/slope_dof)**12}

l2_ylim = [.5e-10,2e-3]
slope_L2_list = {
  "o2": 2e2*(1.0/slope_dof)**2,
  "o3_2": 5e4*(1.0/slope_dof)**3,         
  "o4": 3e7*(1.0/slope_dof)**4,   
  "o7": 1.2e14*(1.0/slope_dof)**7,  
  "o8": 2e16*(1.0/slope_dof)**8,
  "o9": .5e19*(1.0/slope_dof)**9,
  "o12": .2e27*(1.0/slope_dof)**12}

linf_ylim = [1e-8,1e-1]
slope_Linf_list = {
  "o2": 2e4*(1.0/slope_dof)**2,
  "o3": 1.5e7*(1.0/slope_dof)**3,  
  "o4": 4e9*(1.0/slope_dof)**4,   
#  "o4_2": 1e6*(1.0/slope_dof)**4,   
  "o6": 4e13*(1.0/slope_dof)**6,
#  "o8": 2e18*(1.0/slope_dof)**8,
  "o10": 2e24*(1.0/slope_dof)**10,    
  "o12": 5e28*(1.0/slope_dof)**12}


mncs.mkgraph_num_convergence(
  [l1_error_list,l2_error_list,linf_error_list], 
  [slope_dof, slope_dof, slope_dof],  
  [slope_L1_list, slope_L2_list, slope_Linf_list], 
  [f"L1 error", f"L2 error", f"Linf error"], 
  [l1_ylim, l2_ylim, linf_ylim], f"{OUT_DIR}/numerror_{varname}_t2hr.pdf" )

# W --------
varname="W"

l1_error_list = {}; l2_error_list = {}; linf_error_list = {}
for key, numerror in numerror_list.items():
  l1_error_list[key], l2_error_list[key], linf_error_list[key] = mncs.resoldep_numerrordata_sel_time(varname, time, numerror_list[key], dof_list[key])

slope_dof = l1_error_list["P3"].DOF
l1_ylim = [1e-12,5e-4]
slope_L1_list = {
  "o2": .5e2*(1.0/slope_dof)**2,
  "o4": 2e6*(1.0/slope_dof)**4,   
  "o8": 1e16*(1.0/slope_dof)**8,
  "o9_2": .6e18*(1.0/slope_dof)**9,    
  "o12": .3e27*(1.0/slope_dof)**12}

l2_ylim = [3e-11,.5e-2]
slope_L2_list = {
  "o2": 2e2*(1.0/slope_dof)**2,
  "o3": 3e4*(1.0/slope_dof)**3,  
  "o4": 1e7*(1.0/slope_dof)**4,   
  "o7": 1.5e14*(1.0/slope_dof)**7,  
  "o8": .6e17*(1.0/slope_dof)**8,
  "o8_2": 1.3e16*(1.0/slope_dof)**8,  
  "o12": .4e28*(1.0/slope_dof)**12}

linf_ylim = [1e-8,2e-1]
slope_Linf_list = {
  "o2": 3e4*(1.0/slope_dof)**2,
  "o3": 1.2e7*(1.0/slope_dof)**3,  
  "o4": 2e9*(1.0/slope_dof)**4,   
#  "o4_2": 1e6*(1.0/slope_dof)**4,   
  "o6": 5e13*(1.0/slope_dof)**6,
  "o8": 3e18*(1.0/slope_dof)**8,
  "o12": 1e30*(1.0/slope_dof)**12}

mncs.mkgraph_num_convergence(
  [l1_error_list,l2_error_list,linf_error_list], 
  [slope_dof, slope_dof, slope_dof],  
  [slope_L1_list, slope_L2_list, slope_Linf_list], 
  [f"L1 error", f"L2 error", f"Linf error"], 
  [l1_ylim, l2_ylim, linf_ylim], f"{OUT_DIR}/numerror_{varname}_t2hr.pdf" )

# DDENS --------
varname="DDENS"

l1_error_list = {}; l2_error_list = {}; linf_error_list = {}
for key, numerror in numerror_list.items():
  l1_error_list[key], l2_error_list[key], linf_error_list[key] = mncs.resoldep_numerrordata_sel_time(varname, time, numerror_list[key], dof_list[key])

slope_dof = l1_error_list["P3"].DOF
l1_ylim = [2e-14,1e-6]
slope_L1_list = {
  "o2": 1e-1*(1.0/slope_dof)**2,
  "o3_2": 2e1*(1.0/slope_dof)**3,   
  "o4": 1.5e4*(1.0/slope_dof)**4, 
  "o7": .5e11*(1.0/slope_dof)**7,    
  "o8": 1.5e13*(1.0/slope_dof)**8,
  "o10": .6e18*(1.0/slope_dof)**10,    
  "o12": .7e23*(1.0/slope_dof)**12}

l2_ylim = [1e-13,5e-6]
slope_L2_list = {
  "o2": 5e-1*(1.0/slope_dof)**2,
  "o3_2": 9e1*(1.0/slope_dof)**3,   
  "o4": 7e4*(1.0/slope_dof)**4, 
  "o7": 3e11*(1.0/slope_dof)**7,  
  "o8": 1e14*(1.0/slope_dof)**8,
  "o10": 5e18*(1.0/slope_dof)**10,  
  "o12": .5e24*(1.0/slope_dof)**12}

linf_ylim = [3e-11,.5e-3]
slope_Linf_list = {
  "o2": 3e1*(1.0/slope_dof)**2,
  "o3": 1.5e4*(1.0/slope_dof)**3,   
  "o4": 1e7*(1.0/slope_dof)**4, 
  "o6": .7e11*(1.0/slope_dof)**6,  
  "o8": 1e16*(1.0/slope_dof)**8,
  "o9": 2e18*(1.0/slope_dof)**9,  
  "o12": .5e26*(1.0/slope_dof)**12}

mncs.mkgraph_num_convergence(
  [l1_error_list,l2_error_list,linf_error_list], 
  [slope_dof, slope_dof, slope_dof],  
  [slope_L1_list, slope_L2_list, slope_Linf_list], 
  [f"L1 error", f"L2 error", f"Linf error"], 
  [l1_ylim, l2_ylim, linf_ylim], f"{OUT_DIR}/numerror_{varname}_t2hr.pdf" )

# DRHOT --------
varname="THERM"

l1_error_list = {}; l2_error_list = {}; linf_error_list = {}
for key, numerror in numerror_list.items():
  l1_error_list[key], l2_error_list[key], linf_error_list[key] = mncs.resoldep_numerrordata_sel_time(varname, time, numerror_list[key], dof_list[key])

slope_dof = l1_error_list["P3"].DOF
l1_ylim = [1e-12,5e-5]
slope_L1_list = {
  "o2": .5e1*(1.0/slope_dof)**2,
#  "o3": 3.5e0*(1.0/slope_dof)**3,  
  "o4": 4e5*(1.0/slope_dof)**4,   
  "o7": 1e12*(1.0/slope_dof)**7,    
#  "o8": 3e15*(1.0/slope_dof)**8,
  "o10": .8e20*(1.0/slope_dof)**10,
  "o12": .5e26*(1.0/slope_dof)**12}

l2_ylim = [1e-11,3e-4]
slope_L2_list = {
  "o2": 4e1*(1.0/slope_dof)**2,
  "o3": 7e3*(1.0/slope_dof)**3,  
  "o4": 2e6*(1.0/slope_dof)**4,   
  "o7": 3e13*(1.0/slope_dof)**7,    
#  "o8": 2e16*(1.0/slope_dof)**8,
  "o12": 1e27*(1.0/slope_dof)**12}

linf_ylim = [3e-9,5e-2]
slope_Linf_list = {
  "o2": .3e4*(1.0/slope_dof)**2,
  "o3": 1e6*(1.0/slope_dof)**3,      
  "o4": .3e9*(1.0/slope_dof)**4,   
  "o7": 2e16*(1.0/slope_dof)**7,  
#  "o8": 1e19*(1.0/slope_dof)**8,
  "o9": 7e20*(1.0/slope_dof)**9,  
  "o12": 2e29*(1.0/slope_dof)**12}

mncs.mkgraph_num_convergence(
  [l1_error_list,l2_error_list,linf_error_list], 
  [slope_dof, slope_dof, slope_dof],  
  [slope_L1_list, slope_L2_list, slope_Linf_list], 
  [f"L1 error", f"L2 error", f"Linf error"], 
  [l1_ylim, l2_ylim, linf_ylim], f"{OUT_DIR}/numerror_{varname}_t2hr.pdf" )
