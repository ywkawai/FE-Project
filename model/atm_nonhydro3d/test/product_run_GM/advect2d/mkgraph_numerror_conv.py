import os
import mkgraph_numerror_conv_sub as mncs

TOP_DIR="./"
OUT_DIR="./analysis/num_error"

exp_dir_list = {
    "0deg": f"{TOP_DIR}/solid_body_rot_alph_0deg/GAUSSIAN/",
    "45deg": f"{TOP_DIR}/solid_body_rot_alph_45deg/GAUSSIAN/",
    "90deg": f"{TOP_DIR}/solid_body_rot_alph_90deg/GAUSSIAN/"    
}

p1_list = ["Eh16Ez12P1", "Eh32Ez12P1", "Eh64Ez12P1", "Eh128Ez12P1"]
p1_dof_list = {"Eh16Ez12P1":32, "Eh32Ez12P1":64, "Eh64Ez12P1":128, "Eh128Ez12P1":256}

p3_list = ["Eh8Ez6P3", "Eh16Ez6P3", "Eh32Ez6P3", "Eh64Ez6P3"]
p3_dof_list = {"Eh8Ez6P3":32, "Eh16Ez6P3":64, "Eh32Ez6P3":128, "Eh64Ez6P3":256}

p7_list = ["Eh4Ez3P7", "Eh8Ez3P7", "Eh16Ez3P7", "Eh32Ez3P7"]
p7_dof_list = {"Eh4Ez3P7":32, "Eh8Ez3P7":64, "Eh16Ez3P7":128, "Eh32Ez3P7":256}

p11_list = ["Eh2Ez3P11", "Eh4Ez3P11", "Eh8Ez3P11", "Eh16Ez3P11"]#, "Eh32Ez3P11"]#, "Eh8Ez3P7", "Eh16Ez3P7", "Eh32Ez3P7"]
p11_dof_list = {"Eh2Ez3P11":24, "Eh4Ez3P11":48, "Eh8Ez3P11":96, "Eh16Ez3P11":192}#, "Eh32Ez3P11":384}

numerror_list = {}
dof_list = {"P1":p1_dof_list, "P3":p3_dof_list, "P7":p7_dof_list, "P11":p11_dof_list}
for key, exp_dir in exp_dir_list.items():
    numerror_list[f"P1_{key}"] = mncs.get_numerror_data_list(exp_dir, p1_list)
    numerror_list[f"P3_{key}"] = mncs.get_numerror_data_list(exp_dir, p3_list)
    numerror_list[f"P7_{key}"] = mncs.get_numerror_data_list(exp_dir, p7_list)
    numerror_list[f"P11_{key}"] = mncs.get_numerror_data_list(exp_dir, p11_list)

os.makedirs(f"{OUT_DIR}", exist_ok=True)

#-- time = 0 s --------
varname="PTracer"
time=0

l1_error_list = {}; l2_error_list = {}; linf_error_list = {}
for porder in [1, 3, 7, 11]:
  for key, exp_dir in exp_dir_list.items():
    l1_error_list[f"P{porder}_{key}"], l2_error_list[f"P{porder}_{key}"], linf_error_list[f"P{porder}_{key}"] = mncs.resoldep_numerrordata_sel_time(varname, time, numerror_list[f"P{porder}_{key}"], dof_list[f"P{porder}"])

slope_dof = l1_error_list["P7_0deg"].DOF
l1_ylim = [1e-17,1e-3]
slope_L1_list = {
  "o2": 2e-0*(1.0/slope_dof)**2,
  "o4": 2e2*(1.0/slope_dof)**4, 
  "o8": 1e7*(1.0/slope_dof)**8,
  "o12": .6e12*(1.0/slope_dof)**12}

l2_ylim = [1e-16,1e-2]
slope_L2_list = {
  "o2": 1.2e1*(1.0/slope_dof)**2,
  "o4": 5e2*(1.0/slope_dof)**4, 
  "o8": 2e7*(1.0/slope_dof)**8,
  "o12": .8e12*(1.0/slope_dof)**12}

linf_ylim = [1e-14,5e-1]
slope_Linf_list = {
  "o2": 2e2*(1.0/slope_dof)**2,
  "o4": 3e4*(1.0/slope_dof)**4, 
  "o8": 2e9*(1.0/slope_dof)**8,
  "o12": 2e14*(1.0/slope_dof)**12}

mncs.mkgraph_num_convergence(
  [l1_error_list,l2_error_list,linf_error_list], 
  [slope_dof, slope_dof, slope_dof],  
  [slope_L1_list, slope_L2_list, slope_Linf_list], 
  [f"L1 error", f"L2 error", f"Linf error"], 
  [l1_ylim, l2_ylim, linf_ylim], f"{OUT_DIR}/numerror_t0day.png" )

#-- time = 12 * 86400 s --------

varname="PTracer"
time=86400*12

l1_error_list = {}; l2_error_list = {}; linf_error_list = {}
for porder in [1, 3, 7, 11]:
  for key, exp_dir in exp_dir_list.items():
    l1_error_list[f"P{porder}_{key}"], l2_error_list[f"P{porder}_{key}"], linf_error_list[f"P{porder}_{key}"] = mncs.resoldep_numerrordata_sel_time(varname, time, numerror_list[f"P{porder}_{key}"], dof_list[f"P{porder}"])

slope_dof = l1_error_list["P7_0deg"].DOF
l1_ylim = [1e-14,1e-2]
slope_L1_list = {
  "o2": 2e-0*(1.0/slope_dof)**2,
  "o4": 1e2*(1.0/slope_dof)**4, 
  "o8": .7e7*(1.0/slope_dof)**8,
  "o4_2": 1.5e-4*(1.0/slope_dof)**4,  
  "o12": .6e12*(1.0/slope_dof)**12}

l2_ylim = [1e-13,1e-1]
slope_L2_list = {
  "o2": 1.2e1*(1.0/slope_dof)**2,
  "o4": .8e3*(1.0/slope_dof)**4, 
  "o8": 4e7*(1.0/slope_dof)**8,
  "o4_2": 1e-3*(1.0/slope_dof)**4,
  "o12": 3e12*(1.0/slope_dof)**12}

linf_ylim = [1e-12,1e0]
slope_Linf_list = {
  "o2": 2e2*(1.0/slope_dof)**2,
  "o4": 3e4*(1.0/slope_dof)**4, 
  "o8": 2e9*(1.0/slope_dof)**8, 
  "o4_2": 1e-1*(1.0/slope_dof)**4,
  "o12": 1e14*(1.0/slope_dof)**12}


mncs.mkgraph_num_convergence(
  [l1_error_list,l2_error_list,linf_error_list], 
  [slope_dof, slope_dof, slope_dof],  
  [slope_L1_list, slope_L2_list, slope_Linf_list], 
  [f"L1 error", f"L2 error", f"Linf error"], 
  [l1_ylim, l2_ylim, linf_ylim], f"{OUT_DIR}/numerror_t12day.png" )
