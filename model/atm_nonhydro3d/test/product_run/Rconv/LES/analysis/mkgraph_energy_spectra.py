import os
import sys
import numpy as np
sys.path.append(os.path.join(os.path.dirname('__file__'), '../../analysis_lib/'))
import energy_spectra_common as common

#--
OUTDIR="fig_energy_spectra/"
outfigext = "svg"

cr = 0.85

exp_name_list = [
  'Dx25m_P3', 
  'Dx25m_P7', 
  'Dx27m_P11', 
  'Dx12.5m_P3',
  'Dx12.5m_P7', 
  'Dx13m_P11', 
]

dir_ind_list = { 'Dx25m_P3': np.arange(4,25, 1), 'Dx25m_P7': np.arange(4,25, 1), 'Dx27m_P11': np.arange(4,25, 1), 
                 'Dx12.5m_P3': np.arange(4,25, 1), 'Dx12.5m_P7': np.arange(4,25, 1), 'Dx13m_P11': np.arange(4,25, 1), 
                 'Dx3.1m_P7': np.arange(13,101, 1), 
}

K_EXP_NAME_LABEL='Dx3.1m_P7'
YLIM_RANGE=[1e-6,2e0]
SLOPE_m35_ampl=6.5e-5

exp_label_list = {
  "Dx25m_P3":  "Δ=25m,P3",  "Dx25m_P7":  "Δ=25m,P7", "Dx27m_P11":  "Δ=27m,P11",  
  "Dx12.5m_P3":  "Δ=13m,P3", "Dx12.5m_P7":  "Δ=13m,P7", "Dx13m_P11": "Δ=13m,P11", 
  "Dx3.1m_P7": "Δ=3.1m,P7", 
}
exp_color_list = {
  "Dx25m_P3": "black", "Dx25m_P7": "black", "Dx27m_P11":  "black",   
  "Dx12.5m_P3": "red", "Dx12.5m_P7": "red", "Dx13m_P11": "red",
  "Dx3.1m_P7": "blue", 
}
exp_ltype_list = {
  "Dx25m_P3": ":", "Dx25m_P7": "--", "Dx27m_P11": "-", 
  "Dx12.5m_P3": ":", "Dx12.5m_P7": "--", "Dx13m_P11": "-", 
  "Dx3.1m_P7": "--",  
}
exp_ltype_width = {
  "Dx25m_P3": 3, "Dx25m_P7": 3, "Dx27m_P11": 3, 
  "Dx12.5m_P3": 3, "Dx12.5m_P7": 3, "Dx13m_P11": 3, 
  "Dx3.1m_P7": 3,  
}
exp_ncut = {
  "Dx25m_P3": 64, "Dx25m_P7": 64, "Dx27m_P11": 60,  
  "Dx12.5m_P3": 128, "Dx12.5m_P7": 128, "Dx13m_P11": 120, 
  "Dx3.1m_P7": 512, 
}

ZLEVEL_list = [800]

# Reference
exp_name_ref = "Dx3.1m_P7"
tmp_ref_dir = f"../../DNS_3/tmp_data_energy_spectra/tmp_{exp_name_ref}"

#---------------------------------------------------------------------------------------
os.makedirs(OUTDIR, exist_ok=True)

ke_spectra_3dvel_list = {}

for exp_name in exp_name_list:
  tmp_dir = f"../tmp_data_energy_spectra/tmp_{exp_name}"
  print(f"tmp_dir={tmp_dir}")

  common.read_spectra_data( exp_name, tmp_dir, dir_ind_list[exp_name], ke_spectra_3dvel_list, exp_ncut[exp_name] )

# Reference data
common.read_spectra_data( exp_name_ref, tmp_ref_dir, dir_ind_list[exp_name_ref], ke_spectra_3dvel_list, exp_ncut[exp_name_ref] )
exp_name_list_ = exp_name_list.copy()
exp_name_list_.append(exp_name_ref)

for zind, zlev in enumerate(ZLEVEL_list):
  common.create_fig_energy_spectra( ke_spectra_3dvel_list, zlev, f"{OUTDIR}/3dvel_spectra_z{zlev}m.{outfigext}", 
                                   exp_ltype_list, exp_ltype_width, exp_color_list, exp_label_list, YLIM_RANGE, K_EXP_NAME_LABEL,  
                                   SLOPE_m35_ampl ) # 1.7e-5  
  common.create_fig_energy_spectra_diffm53( ke_spectra_3dvel_list, zlev, f"{OUTDIR}/3dvel_spectra_z{zlev}m_diff.{outfigext}", 
                                    exp_ltype_list, exp_ltype_width, exp_color_list, exp_label_list, 
                                    [2e-1,2e0], SLOPE_m35_ampl )  

  common.create_fig_energy_spectra_diff_refexp( ke_spectra_3dvel_list, K_EXP_NAME_LABEL, zlev, cr, 
                                             exp_ncut, exp_name_list_, exp_label_list, exp_color_list, exp_ltype_list, exp_ltype_width,
                                             f"{OUTDIR}/3dvel_spectra_diff_refexp_z{ZLEVEL_list[0]}m.{outfigext}" )