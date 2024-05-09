import os
import mkgraph_numerror_conv_sub as mncs

TOP_DIR="./"
OUT_DIR="./analysis_out/num_error"

p1_list = ["Eh16Ez12P1", "Eh32Ez24P1", "Eh64Ez48P1", "Eh128Ez96P1"]
p1_dof_list = {"Eh16Ez12P1":32, "Eh32Ez24P1":64, "Eh64Ez48P1":128, "Eh128Ez96P1":256}

p3_list = ["Eh8Ez6P3", "Eh16Ez12P3", "Eh32Ez24P3", "Eh64Ez48P3"]
p3_dof_list = {"Eh8Ez6P3":32, "Eh16Ez12P3":64, "Eh32Ez24P3": 128, "Eh64Ez48P3": 256}

p7_list = ["Eh4Ez3P7", "Eh8Ez6P7", "Eh16Ez12P7", "Eh32Ez24P7"]
p7_dof_list = {"Eh4Ez3P7":32, "Eh8Ez6P7":64, "Eh16Ez12P7":128, "Eh32Ez24P7":256}

p7_dtx0p5_list = ["Eh4Ez3P7_dtx0.5", "Eh8Ez6P7_dtx0.5", "Eh16Ez12P7_dtx0.5", "Eh32Ez24P7_dtx0.5"]
p7_dtx0p5_dof_list = {"Eh4Ez3P7_dtx0.5":32, "Eh8Ez6P7_dtx0.5":64, "Eh16Ez12P7_dtx0.5":128, "Eh32Ez24P7_dtx0.5":256}

p7_dtx0p25_list = ["Eh4Ez3P7_dtx0.25", "Eh8Ez6P7_dtx0.25", "Eh16Ez12P7_dtx0.25", "Eh32Ez24P7_dtx0.25"]
p7_dtx0p25_dof_list = {"Eh4Ez3P7_dtx0.25":32, "Eh8Ez6P7_dtx0.25":64, "Eh16Ez12P7_dtx0.25":128, "Eh32Ez24P7_dtx0.25":256}

p11_list = ["Eh3Ez2P11", "Eh6Ez4P11", "Eh12Ez8P11"]
p11_dof_list = {"Eh3Ez2P11":36, "Eh6Ez4P11":72, "Eh12Ez8P11":144}

p11_dtx0p5_list = ["Eh3Ez2P11_dtx0.5", "Eh6Ez4P11_dtx0.5", "Eh12Ez8P11_dtx0.5"]
p11_dtx0p5_dof_list = {"Eh3Ez2P11_dtx0.5":36, "Eh6Ez4P11_dtx0.5":72, "Eh12Ez8P11_dtx0.5":144}

p11_dtx0p25_list = ["Eh3Ez2P11_dtx0.25", "Eh6Ez4P11_dtx0.25", "Eh12Ez8P11_dtx0.25"]
p11_dtx0p25_dof_list = {"Eh3Ez2P11_dtx0.25":36, "Eh6Ez4P11_dtx0.25":72, "Eh12Ez8P11_dtx0.25":144}

#------

numerror_list = {}
dof_list = {"P1":p1_dof_list, "P3":p3_dof_list, 
            "P7":p7_dof_list, "P7_dtx0.5":p7_dtx0p5_dof_list, "P7_dtx0.25":p7_dtx0p25_dof_list, 
            "P11":p11_dof_list, "P11_dtx0.5":p11_dtx0p5_dof_list, "P11_dtx0.25":p11_dtx0p25_dof_list

          }
exp_dir="rhot_hevi"
analysis_data="analysis_3/NUMERROR_LOG.peall"
numerror_list[f"P1"] = mncs.get_numerror_data_list(exp_dir, p1_list, analysis_data)
numerror_list[f"P3"] = mncs.get_numerror_data_list(exp_dir, p3_list, analysis_data)
numerror_list[f"P7"] = mncs.get_numerror_data_list(exp_dir, p7_list, analysis_data)
numerror_list[f"P7_dtx0.5"] = mncs.get_numerror_data_list(exp_dir, p7_dtx0p5_list, analysis_data)
numerror_list[f"P7_dtx0.25"] = mncs.get_numerror_data_list(exp_dir, p7_dtx0p25_list, analysis_data)
numerror_list[f"P11"] = mncs.get_numerror_data_list(exp_dir, p11_list, analysis_data)
numerror_list[f"P11_dtx0.5"] = mncs.get_numerror_data_list(exp_dir, p11_dtx0p5_list, analysis_data)
numerror_list[f"P11_dtx0.25"] = mncs.get_numerror_data_list(exp_dir, p11_dtx0p25_list, analysis_data)

os.makedirs(f"{OUT_DIR}", exist_ok=True)


#-- time = 2 * 86400 s 
time=172800

# DENS --------
varname="DDENS"

l1_error_list = {}; l2_error_list = {}; linf_error_list = {}
for key, numerror in numerror_list.items():
  l1_error_list[key], l2_error_list[key], linf_error_list[key] = mncs.resoldep_numerrordata_sel_time(varname, time, numerror_list[key], dof_list[key])

slope_dof = l1_error_list["P7"].DOF
l1_ylim = [4e-14,1e-6]
slope_L1_list = {
  "o2": 1e-3*(1.0/slope_dof)**2,
  "o4": 2.5e-2*(1.0/slope_dof)**4, 
  "o3_2": 6e-6*(1.0/slope_dof)**3,       
  "o8": 8e2*(1.0/slope_dof)**8,
  "o12": 1e8*(1.0/slope_dof)**12}

l2_ylim = [1e-13,2e-6]
slope_L2_list = {
  "o2": 2e-3*(1.0/slope_dof)**2,
  "o4": 4e-2*(1.0/slope_dof)**4, 
  "o3_2": 1e-5*(1.0/slope_dof)**3,     
  "o8": 1.5e3*(1.0/slope_dof)**8,
  "o12": 5e8*(1.0/slope_dof)**12}

linf_ylim = [3e-11,5e-6]
slope_Linf_list = {
  "o1": .3e-3*(1.0/slope_dof)**1,  
  "o2": 1e-2*(1.0/slope_dof)**2,
  "o3_1": 6e-3*(1.0/slope_dof)**3,   
  "o4": 3e-1*(1.0/slope_dof)**4, 
#  "o3_2": 5e-4*(1.0/slope_dof)**3,     
  "o8": 2e4*(1.0/slope_dof)**8,
  "o12": 1e10*(1.0/slope_dof)**12}

mncs.mkgraph_num_convergence(
  [l1_error_list,l2_error_list,linf_error_list], 
  [slope_dof, slope_dof, slope_dof],  
  [slope_L1_list, slope_L2_list, slope_Linf_list], 
  [f"L1 error", f"L2 error", f"Linf error"], 
  [l1_ylim, l2_ylim, linf_ylim], f"{OUT_DIR}/numerror_{varname}_t2day.png" )

# U --------
varname="U"

l1_error_list = {}; l2_error_list = {}; linf_error_list = {}
for key, numerror in numerror_list.items():
  l1_error_list[key], l2_error_list[key], linf_error_list[key] = mncs.resoldep_numerrordata_sel_time(varname, time, numerror_list[key], dof_list[key])

slope_dof = l1_error_list["P7"].DOF
l1_ylim = [1e-18,7e-11]
slope_L1_list = {
  "o2": 7e-8*(1.0/slope_dof)**2,
  "o4": 2e-6*(1.0/slope_dof)**4, 
  "o3_2": 5e-11*(1.0/slope_dof)**3,    
  "o8": 6e-2*(1.0/slope_dof)**8,
  "o12": 5e3*(1.0/slope_dof)**12}

l2_ylim = [2e-18,1.2e-10]
slope_L2_list = {
  "o2": 2e-7*(1.0/slope_dof)**2,
  "o4": 4e-6*(1.0/slope_dof)**4, 
  "o3_2": 1.1e-10*(1.0/slope_dof)**3,      
  "o8": 1.2e-1*(1.0/slope_dof)**8,
  "o12": 1e4*(1.0/slope_dof)**12
  }

linf_ylim = [1.5e-15,2e-9]
slope_Linf_list = {
  "o2": 2e-6*(1.0/slope_dof)**2,
  "o4": 5e-5*(1.0/slope_dof)**4, 
  "o8": 3e-0*(1.0/slope_dof)**8,
  "o12": 1e5*(1.0/slope_dof)**12}

mncs.mkgraph_num_convergence(
  [l1_error_list,l2_error_list,linf_error_list], 
  [slope_dof, slope_dof, slope_dof],  
  [slope_L1_list, slope_L2_list, slope_Linf_list], 
  [f"L1 error", f"L2 error", f"Linf error"], 
  [l1_ylim, l2_ylim, linf_ylim], f"{OUT_DIR}/numerror_{varname}_t2day.png" )

# W --------
varname="W"

l1_error_list = {}; l2_error_list = {}; linf_error_list = {}
for key, numerror in numerror_list.items():
  l1_error_list[key], l2_error_list[key], linf_error_list[key] = mncs.resoldep_numerrordata_sel_time(varname, time, numerror_list[key], dof_list[key])

slope_dof = l1_error_list["P7"].DOF
l1_ylim = [1e-11,4e-6]
slope_L1_list = {
  "o2": 5e-3*(1.0/slope_dof)**2,  
  "o4": 3.5e-1*(1.0/slope_dof)**4, 
  "o2_2": 1.6e-5*(1.0/slope_dof)**2,  
  "o3_2": 7e-5*(1.0/slope_dof)**3,   
  "o8": 2e4*(1.0/slope_dof)**8,
  "o12": 1e9*(1.0/slope_dof)**12}

l2_ylim = [2e-11,5e-6]
slope_L2_list = {
  "o1": 1.2e-4*(1.0/slope_dof)**1,    
  "o2": 9e-3*(1.0/slope_dof)**2,  
  "o4": 8e-1*(1.0/slope_dof)**4, 
  "o2_2": 5e-5*(1.0/slope_dof)**2,  
  "o3_2": 1.2e-4*(1.0/slope_dof)**3,    
  "o8": 5e4*(1.0/slope_dof)**8,
#  "o12": 1e9*(1.0/slope_dof)**12
  }

linf_ylim = [1e-10,3e-5]
slope_Linf_list = {
  "o1": 7e-4*(1.0/slope_dof)**1,  
  "o2": 8e-2*(1.0/slope_dof)**2,
  "o3_2": 1.8e-1*(1.0/slope_dof)**3,   
  "o4": 9e-0*(1.0/slope_dof)**4, 
  "o2_2": 1.5e-4*(1.0/slope_dof)**2, 
  "o3_2": 1.5e-3*(1.0/slope_dof)**3,         
  "o8": 1e6*(1.0/slope_dof)**8,
#  "o12": 1.4e11*(1.0/slope_dof)**12
  }

mncs.mkgraph_num_convergence(
  [l1_error_list,l2_error_list,linf_error_list], 
  [slope_dof, slope_dof, slope_dof],  
  [slope_L1_list, slope_L2_list, slope_Linf_list], 
  [f"L1 error", f"L2 error", f"Linf error"], 
  [l1_ylim, l2_ylim, linf_ylim], f"{OUT_DIR}/numerror_{varname}_t2day.png" )

# W --------
varname="DRHOT"

l1_error_list = {}; l2_error_list = {}; linf_error_list = {}
for key, numerror in numerror_list.items():
  l1_error_list[key], l2_error_list[key], linf_error_list[key] = mncs.resoldep_numerrordata_sel_time(varname, time, numerror_list[key], dof_list[key])

slope_dof = l1_error_list["P7"].DOF
l1_ylim = [5e-12,6e-5]
slope_L1_list = {
  "o2": 8e-2*(1.0/slope_dof)**2,
  "o4": 2e-0*(1.0/slope_dof)**4, 
  "o3_2": 2e-3*(1.0/slope_dof)**3,   
  "o8": 8e4*(1.0/slope_dof)**8,
  "o12": 1e10*(1.0/slope_dof)**12}

l2_ylim = [5e-12,1.2e-4]
slope_L2_list = {
  "o2": 1.5e-1*(1.0/slope_dof)**2,
  "o4": 3e-0*(1.0/slope_dof)**4, 
  "o3_2": 4e-3*(1.0/slope_dof)**3,     
  "o8": 1.1e5*(1.0/slope_dof)**8,
  "o12": 2e10*(1.0/slope_dof)**12}

linf_ylim = [6e-12,6e-4]
slope_Linf_list = {
  "o2": 8e-1*(1.0/slope_dof)**2,
  "o4": 3e1*(1.0/slope_dof)**4, 
  "o3_2": 1.5e-2*(1.0/slope_dof)**3,  
  "o8": 2e6*(1.0/slope_dof)**8,
  "o12": 3e11*(1.0/slope_dof)**12}

mncs.mkgraph_num_convergence(
  [l1_error_list,l2_error_list,linf_error_list], 
  [slope_dof, slope_dof, slope_dof],  
  [slope_L1_list, slope_L2_list, slope_Linf_list], 
  [f"L1 error", f"L2 error", f"Linf error"], 
  [l1_ylim, l2_ylim, linf_ylim], f"{OUT_DIR}/numerror_{varname}_t2day.png" )