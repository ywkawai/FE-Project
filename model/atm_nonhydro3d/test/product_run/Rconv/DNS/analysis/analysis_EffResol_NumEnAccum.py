import os
import sys
import numpy as np
sys.path.append(os.path.join(os.path.dirname('__file__'), '../../analysis_lib/'))
import energy_spectra_common as common
import eval_EffResol_NumEnAccum as analysis
#--

OUTDIR="comp_cost_info/"

ref_exp_name ="Dx3.1m_P7"
cr = 0.85
cr2 = 1.05

exp_name_list = [
  'Dx25m_P3', 
  'Dx25m_P7', 
  'Dx27m_P11', 
  'Dx12.5m_P3',
  'Dx12.5m_P7', 
  'Dx13m_P11', 
  'Dx6.3m_P3', 
  'Dx6.3m_P7', 
  'Dx6.7m_P11',
  'Dx3.1m_P7',           
]

dir_ind_list = { 'Dx25m_P3': np.arange(4,25, 1), 'Dx25m_P7': np.arange(4,25, 1), 'Dx27m_P11': np.arange(4,25, 1), 
                 'Dx12.5m_P3': np.arange(8,29, 1), 'Dx12.5m_P7': np.arange(4,25, 1), 'Dx13m_P11': np.arange(8,29, 1), 
                 'Dx6.3m_P3': np.arange(13,58, 1), 'Dx6.3m_P7': np.arange(8,53, 1), 'Dx6.7m_P11': np.arange(13,35, 1),
                 'Dx3.1m_P7': np.arange(13,101, 1), 'Dx3.1m_P7_old': np.arange(21,90, 1), 
}
nx_list = {'Dx25m_P3': 128, 'Dx25m_P3_MF': 128, 
           'Dx25m_P7': 128, 'Dx25m_P7_2': 128,  'Dx25m_P7_3': 128,  
           'Dx27m_P11': 120, 'Dx27m_P11_MF': 120, 
           'Dx12.5m_P3':256, 'Dx12.5m_P3_MFweak':256, 
           'Dx12.5m_P7': 256, 'Dx12.5m_P7_MF': 256, 
           'Dx13m_P11': 240, 
           'Dx6.3m_P3': 512, 'Dx6.3m_P7': 512, 'Dx6.7m_P11': 480, 
           'Dx3.1m_P7': 1024 }
exp_ncut = {
  "Dx25m_P3": 64, "Dx25m_P3_MF": 64, 
  "Dx25m_P7": 64, "Dx25m_P7_2": 64, "Dx25m_P7_3": 64, 
  "Dx27m_P11": 60, "Dx27m_P11_MF": 60, 
  "Dx12.5m_P3": 128, 'Dx12.5m_P3_MFweak': 128, 
  "Dx12.5m_P7": 128, "Dx12.5m_P7_MF": 128, 
  'Dx13m_P11': 120, 
  'Dx6.3m_P3':256, 'Dx6.3m_P7': 256, 'Dx6.7m_P11': 240, 
  'Dx3.1m_P7': 512, 
}


ZLEVEL_list = [800]
Lx = 3.2e3

#---------------------------------------------------------------------------------------
ke_spectra_3dvel_list = {}
ke_spectra_ratio_list = {}

eff_resol_list = {}
energy_accum_list = {}

#--
for exp_name in exp_name_list:
    tmp_dir = f"../tmp_data_energy_spectra/tmp_{exp_name}"
    print(f"tmp_dir={tmp_dir}")

    common.read_spectra_data( exp_name, tmp_dir, dir_ind_list[exp_name], ke_spectra_3dvel_list, exp_ncut[exp_name] )

#--
ke_spectra_ref = ke_spectra_3dvel_list[ref_exp_name].isel(z=0)

for exp_name in exp_name_list:        
    ncut = np.min( [exp_ncut[exp_name], exp_ncut[ref_exp_name]] ) - 1
    ke_spectra =  ke_spectra_3dvel_list[exp_name].isel(z=0)
    k = ke_spectra.k[0:ncut]
    ke_spectra_ratio = ke_spectra[0:ncut]/ke_spectra_ref[0:ncut]
    ke_spectra_ratio_list[exp_name] = ke_spectra_ratio

#-
for exp_name in exp_name_list:
    if exp_name == ref_exp_name:
        break
    ke_spectra_ratio =  ke_spectra_ratio_list[exp_name]
    k_eff, eff_resol_list[exp_name] = analysis.eval_effective_resolution( ke_spectra_ratio[5:-5], 
                                                                         cr, Lx/float(nx_list[exp_name]) )
    
    k1 = 2.0*np.pi/800.0; k2 = k_eff
    # k1 = 2.0*np.pi/1600.0; k2 = k_eff
    ke_spectra_ = ke_spectra_3dvel_list[exp_name].isel(z=0)
    abs_error, rel_error = analysis.eval_energy_pile_error(ke_spectra_, ke_spectra_ref, k1, k2, cr=cr2)
    energy_accum_list[exp_name] = rel_error * 100.0
     
    print(f"{exp_name} : {2.0*np.pi/k_eff:12.1f} {eff_resol_list[exp_name]:12.1f}  {100*rel_error:12.1f} %")

#-
num_energy_accum_equiv_dx_list = analysis.eval_energy_pile_error_equiv_dx(energy_accum_list)

os.makedirs(OUTDIR, exist_ok=True)
analysis.output_energy_spectra_eval(eff_resol_list, energy_accum_list, num_energy_accum_equiv_dx_list, 
                                    f"{OUTDIR}/EffResol_NumEnAccum.dat")
