import os
import sys
import numpy as np
sys.path.append(os.path.join(os.path.dirname('__file__'), '../../analysis_lib/'))
import energy_spectra_common as common

#--
exp_top_dir="../"

exp_name_list = [
  # 'Dx25m_P3',   # Start: 4-24
  # 'Dx25m_P7',      # Start: 4-24
  'Dx27m_P11',     # Start: 4-24
  # 'Dx12.5m_P3', # Start: 8-24
  #  'Dx12.5m_P7',   # Start: 4-24
  #  'Dx13m_P11'      # Start: 8-24
  # 'Dx6.3m_P3',      # Start: 13-57
  # 'Dx6.3m_P7',      # Start: 8-52
  # 'Dx6.7m_P11',      # Start: 13-34
  # 'Dx3.1m_P7'      # Start: 13-100
  ]
dir_ind_list = { 
               'Dx25m_P3': np.arange(21,25, 1), 'Dx25m_P7': np.arange(22,25, 1), 'Dx27m_P11': np.arange(24,25, 1),  
                'Dx12.5m_P3': np.arange(23,25, 1), 'Dx12.5m_P7': np.arange(24,25, 1), 'Dx13m_P11': np.arange(22,25, 1), 
                'Dx6.3m_P3': np.arange(55,58, 1), 'Dx6.3m_P7': np.arange(49,53, 1), 'Dx6.7m_P11': np.arange(32,35, 1),
                'Dx3.1m_P7': np.arange(91,101, 1) 
              }
runno_inidata_list = { 'Dx25m_P3': 2, 'Dx25m_P7': 1, 'Dx27m_P11': 2, 
                 'Dx12.5m_P3': 7, 'Dx12.5m_P7': 1, 'Dx13m_P11': 7, 
                 'Dx6.3m_P3': 12, 'Dx6.3m_P7': 7, 'Dx6.7m_P11': 12,
                 'Dx3.1m_P7': 12
              }
time_list_eachrun_list = { 
                'Dx25m_P3': np.arange(0, 7200, 900), 'Dx25m_P7': np.arange(0, 7200, 900), 'Dx27m_P11': np.arange(0, 7200, 900), 
                'Dx12.5m_P3': np.arange(0, 7200, 900), 'Dx12.5m_P7': np.arange(0, 7200, 900), 'Dx13m_P11': np.arange(0, 7200, 900), 
                'Dx6.3m_P3': np.arange(0, 3600, 900), 'Dx6.3m_P7': np.arange(0, 3600, 900), 'Dx6.7m_P11': np.arange(0, 7200, 900),
                'Dx3.1m_P7': np.arange(0, 1800, 900)
}
penum_list = {'Dx25m_P3': 64, 'Dx25m_P7': 64, 'Dx27m_P11': 100, 
              'Dx12.5m_P3': 64, 'Dx12.5m_P7': 64, 'Dx13m_P11': 100,
              'Dx6.3m_P3': 64, 'Dx6.3m_P7': 64, 'Dx6.7m_P11': 400, 
              'Dx3.1m_P7': 128, 
              }
nx_list = {'Dx25m_P3': 128, 'Dx25m_P7': 128, 'Dx27m_P11': 120,
           'Dx12.5m_P3':256, 'Dx12.5m_P7': 256, 'Dx13m_P11': 240, 
           'Dx6.3m_P3': 512, 'Dx6.3m_P7': 512, 'Dx6.7m_P11': 480, 
           'Dx3.1m_P7': 1024 }

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