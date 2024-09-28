import os
import mkgraph_numerror_conv_sub as mncs

TOP_DIR="./"
EXP_DIR_POSTFIX=""
OUT_DIR=f"./analysis/num_error/slotted_cylinders{EXP_DIR_POSTFIX}"

exp_dir_list = {
    "Cyl": f"{TOP_DIR}/slotted_cylinders{EXP_DIR_POSTFIX}/",
}

p1_list = ["Eh16Ez12P1", "Eh32Ez12P1", "Eh64Ez12P1", "Eh128Ez12P1"]
p1_dof_list = {"Eh16Ez12P1":32, "Eh32Ez12P1":64, "Eh64Ez12P1":128, "Eh128Ez12P1":256}

p3_list = ["Eh8Ez6P3", "Eh16Ez6P3", "Eh32Ez6P3", "Eh64Ez6P3"]
p3_dof_list = {"Eh8Ez6P3":32, "Eh16Ez6P3":64, "Eh32Ez6P3":128, "Eh64Ez6P3":256}

p7_list = ["Eh4Ez3P7", "Eh8Ez3P7", "Eh16Ez3P7", "Eh32Ez3P7"]
p7_dof_list = {"Eh4Ez3P7":32, "Eh8Ez3P7":64, "Eh16Ez3P7":128, "Eh32Ez3P7":256}

p7_MF_list = ["Eh4Ez3P7_MF", "Eh8Ez3P7_MF", "Eh16Ez3P7_MF", "Eh32Ez3P7_MF"]
p7_MF_dof_list = {"Eh4Ez3P7_MF":32, "Eh8Ez3P7_MF":64, "Eh16Ez3P7_MF":128, "Eh32Ez3P7_MF":256}

p11_list = ["Eh2Ez3P11", "Eh4Ez3P11", "Eh8Ez3P11", "Eh16Ez3P11"]#, "Eh32Ez3P11"]#, "Eh8Ez3P7", "Eh16Ez3P7", "Eh32Ez3P7"]
p11_dof_list = {"Eh2Ez3P11":24, "Eh4Ez3P11":48, "Eh8Ez3P11":96, "Eh16Ez3P11":192}#, "Eh32Ez3P11":384}

p11_MF_list = ["Eh2Ez3P11_MF", "Eh4Ez3P11_MF", "Eh8Ez3P11_MF", "Eh16Ez3P11_MF"]#, "Eh32Ez3P11"]#, "Eh8Ez3P7", "Eh16Ez3P7", "Eh32Ez3P7"]
p11_MF_dof_list = {"Eh2Ez3P11_MF":24, "Eh4Ez3P11_MF":48, "Eh8Ez3P11_MF":96, "Eh16Ez3P11_MF":192}#, "Eh32Ez3P11":384}

numerror_list = {}
dof_list = {"P1":p1_dof_list, "P3":p3_dof_list, "P7":p7_dof_list, "P11":p11_dof_list}
for key, exp_dir in exp_dir_list.items():
    numerror_list[f"P1_{key}"] = mncs.get_numerror_data_list(exp_dir, p1_list)
    numerror_list[f"P3_{key}"] = mncs.get_numerror_data_list(exp_dir, p3_list)
    numerror_list[f"P7_{key}"] = mncs.get_numerror_data_list(exp_dir, p7_list)
    numerror_list[f"P7_{key}_MF"] = mncs.get_numerror_data_list(exp_dir, p7_MF_list)    
    numerror_list[f"P11_{key}"] = mncs.get_numerror_data_list(exp_dir, p11_list)
    numerror_list[f"P11_{key}_MF"] = mncs.get_numerror_data_list(exp_dir, p11_MF_list)    

os.makedirs(f"{OUT_DIR}", exist_ok=True)

#-- time = 12 * 86400 s --------

varname="PTracer"
time=86400*12

l1_error_list = {}; l2_error_list = {}; linf_error_list = {}
for porder in [1, 3, 7, 11]:
  for key, exp_dir in exp_dir_list.items():
    l1_error_list[f"P{porder}_{key}"], l2_error_list[f"P{porder}_{key}"], linf_error_list[f"P{porder}_{key}"] = mncs.resoldep_numerrordata_sel_time(varname, time, numerror_list[f"P{porder}_{key}"], dof_list[f"P{porder}"])
for porder in [7, 11]:
  for key, exp_dir in exp_dir_list.items():
    l1_error_list[f"P{porder}_{key}_MF"], l2_error_list[f"P{porder}_{key}_MF"], linf_error_list[f"P{porder}_{key}_MF"] = mncs.resoldep_numerrordata_sel_time(varname, time, numerror_list[f"P{porder}_{key}_MF"], dof_list[f"P{porder}"])

for key in l1_error_list.keys():
  # normalize
  l1_error_list[key] = l1_error_list[key] / 1.024e-1
  l2_error_list[key] = l2_error_list[key] / 3.120e-1

slope_dof = l1_error_list["P7_Cyl"].DOF
l1_ylim = [2e-2,1.2e0]
slope_L1_list = {
  "o0.5": 5e0*(1.0/slope_dof)**0.5,    
  "o1": 2e1*(1.0/slope_dof)**1,  
  "o2": 5e2*(1.0/slope_dof)**2,
}

l2_ylim = [2e-2,1.2e0]
slope_L2_list = {
  "o0.5": 4.5e0*(1.0/slope_dof)**0.5,    
  "o0.7": 4.5e0*(1.0/slope_dof)**0.7,      
  "o1": 1.1e1*(1.0/slope_dof)**1,  
  "o2": 3e2*(1.0/slope_dof)**2,
}

linf_ylim = [2e-2,.5e1]
slope_Linf_list = {
  "o0.2": 3e0*(1.0/slope_dof)**0.2,     
#  "o0.3": 3.7e0*(1.0/slope_dof)**0.3,    
  "o0.5": 5e0*(1.0/slope_dof)**0.5,    
  "o1": 2.8e1*(1.0/slope_dof)**1,  
  "o2": 7e2*(1.0/slope_dof)**2,
}

mncs.mkgraph_num_convergence(
  [l1_error_list,l2_error_list,linf_error_list], 
  [slope_dof, slope_dof, slope_dof],  
  [slope_L1_list, slope_L2_list, slope_Linf_list], 
  [f"L1 error", f"L2 error", f"Linf error"], 
#  [l1_ylim, l2_ylim, linf_ylim], f"{OUT_DIR}/numerror_t12day.png" )
  [l1_ylim, l2_ylim, linf_ylim], f"{OUT_DIR}/numerror_t12day.pdf" )
