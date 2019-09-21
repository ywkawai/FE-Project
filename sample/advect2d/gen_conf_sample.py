#!/usr/bin/env python

# Append the path of advect2d test directory as follows. 
#  sys.path.append("FE_PROJECT_TOPDIR/sample/advect2d")
from advect2d_expcase import exp_prepair

TOP_DIR="./"
Ne_List = [4, 8, 16, 32, 64, 128]
PolyOrder_List = [1, 2, 3, 4, 5]
DELT_List = [ [0.0005]*len(Ne_List) ]*len(PolyOrder_List) 
iniGalerkinProjFlag = False

ORG_PE_PATH="/data5/ykawai/FE-Project/sample/advect2d/test_advect2d"

vconst_params = [1.0, 1.0, 0.0, 0.0]
vrigid_body_rot_params = [0.5, 0.5, 1.0, 0.0]
vswirling_params = [0.0, 0.0, 5.0, 0.0]

expcase_params_dict = {
  'Isin_Vcosntant': ['sin', [1.0, 1.0, 0.0, 0.0],  'constant', vconst_params ], 
  'Igaussian-hill_Vcosntant': ['gaussian-hill', [0.5, 0.5, 0.025, 0.0],  'constant', vconst_params ], 
  'Icosine-bell_Vcosntant': ['cosine-bell', [0.5, 0.5, 0.05, 0.0],  'constant', vconst_params ], 
  'Itop-hat_Vcosntant': ['top-hat', [0.5, 0.5, 0.05, 0.0],  'constant', vconst_params ], 
#---
  'Igaussian-hill_Vrigid-body-rot': ['gaussian-hill', [0.25, 0.25, 0.025, 0.0],  'rigid-body-rot', vrigid_body_rot_params ], 
  'Icosine-bell_Vrigid-body-rot': ['cosine-bell', [0.25, 0.25, 0.05, 0.0],  'rigid-body-rot', vrigid_body_rot_params ], 
  'Itop-hat_Vrigid-body-rot': ['top-hat', [0.25, 0.25, 0.05, 0.0],  'rigid-body-rot', vrigid_body_rot_params ], 
#---
  'Igaussian-hill_Vswirling': ['gaussian-hill', [0.25, 0.25, 0.025, 0.0],  'swirling', vswirling_params ], 
  'Icosine-bell_Vswirling': ['cosine-bell', [0.25, 0.25, 0.05, 0.0],  'swirling', vswirling_params ]
}

exp_prepair(TOP_DIR, Ne_List, PolyOrder_List, DELT_List, ORG_PE_PATH, expcase_params_dict, iniGalerkinProjFlag)
