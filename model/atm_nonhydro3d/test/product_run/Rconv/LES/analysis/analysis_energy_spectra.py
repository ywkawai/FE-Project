import os
import sys
import numpy as np
sys.path.append(os.path.join(os.path.dirname('__file__'), '../../analysis_lib/'))
import energy_spectra_common as common

#--
exp_top_dir="../"

exp_name_list = [
  'Dx25m_P3',   # Start: 4-20
  # 'Dx25m_P7',      # Start: 4-20
  # 'Dx27m_P11',     # Start: 4-20
  # 'Dx12.5m_P3', # Start: 4-20
  # 'Dx12.5m_P7',   # Start: 4-20
  # 'Dx13m_P11'      # Start: 4-20
  ]
dir_ind_list = { 
                'Dx25m_P3': np.arange(4,21, 1), 'Dx25m_P7': np.arange(4,21, 1), 'Dx27m_P11': np.arange(4,21, 1), 
                'Dx12.5m_P3': np.arange(13,21, 1), 'Dx12.5m_P7': np.arange(12,21, 1),  'Dx13m_P11': np.arange(13,21, 1), 
              }
runno_inidata_list = { 'Dx25m_P3': 2, 'Dx25m_P7': 1, 'Dx27m_P11': 2, 
                 'Dx12.5m_P3': 2, 'Dx12.5m_P7': 1, 'Dx13m_P11': 2, 
              }
time_list_eachrun_list = { 
                'Dx25m_P3': np.arange(0, 7200, 900), 'Dx25m_P7': np.arange(0, 7200, 900), 'Dx27m_P11': np.arange(0, 7200, 900), 
                'Dx12.5m_P3': np.arange(0, 7200, 900), 'Dx12.5m_P7': np.arange(0, 7200, 900), 'Dx13m_P11': np.arange(0, 7200, 900), 
}
penum_list = {'Dx25m_P3': 64, 'Dx25m_P7': 64, 'Dx27m_P11': 100, 
              'Dx12.5m_P3': 64, 'Dx12.5m_P7': 64, 'Dx13m_P11': 100,
              }
nx_list = {'Dx25m_P3': 128, 'Dx25m_P7': 128, 'Dx27m_P11': 120,
           'Dx12.5m_P3':256, 'Dx12.5m_P7': 256, 'Dx13m_P11': 240, }

ZLEVEL_list = [800]
L = 3.2e3

nproc=8; gem_tmp_data_skip_flag = False
#------------------------------

for exp_name in exp_name_list:
    TIME_list={}
    for runno in dir_ind_list[exp_name]:
        TIME_list[runno] = time_list_eachrun_list[exp_name]    
    common.energy_spectra_analysis(exp_top_dir, exp_name, nx_list[exp_name], penum_list[exp_name], L, 
                                dir_ind_list[exp_name], runno_inidata_list[exp_name], 
                                TIME_list, ZLEVEL_list, 
                                "", nproc, gem_tmp_data_skip_flag)