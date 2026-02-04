import os
import sys
sys.path.append(os.path.join(os.path.dirname('__file__'), '../../analysis_lib/'))
import energy_tseries_common as common
#----

DATA_DIR=".."; ds_list = {}
dx25m_P7 = common.get_moni_data_list(f"{DATA_DIR}/Dx25m_P7/", 24, 3.125e-1)
ds_list["dx25m_P3"] = common.get_moni_data_list(f"{DATA_DIR}/Dx25m_P3/", 22, 5e-1, 3, 2*7200, dx25m_P7.sel(sec=2*7200))
ds_list["dx25m_P7"] = dx25m_P7
ds_list["dx27m_P11"] = common.get_moni_data_list(f"{DATA_DIR}/Dx27m_P11/", 22, 2e-1, 3, 2*7200, dx25m_P7.sel(sec=2*7200))

dx12p5m_P7 = common.get_moni_data_list(f"{DATA_DIR}/Dx12.5m_P7/", 24, 1.5625e-1)
ds_list["dx12.5m_P3"] = common.get_moni_data_list(f"{DATA_DIR}/Dx12.5m_P3/", 17, 2.5e-1, 8, 7*7200, dx12p5m_P7.sel(sec=7*7200))
ds_list["dx12.5m_P7"] = dx12p5m_P7
ds_list["dx13m_P11"] = common.get_moni_data_list(f"{DATA_DIR}/Dx13m_P11/", 17, 1.25e-1, 8, 7*7200, dx12p5m_P7.sel(sec=7*7200))

dx6p3m_P7 = common.get_moni_data_list(f"{DATA_DIR}/Dx6.3m_P7/", 45, 7.8125e-2, 8, 7*7200, ds_list["dx12.5m_P7"].sel(sec=7*7200))
ds_list["dx6.3m_P3"] = common.get_moni_data_list(f"{DATA_DIR}/Dx6.3m_P3/", 42, 1.25e-1, 13, 19*3600, dx6p3m_P7.sel(sec=19*3600))
ds_list["dx6.3m_P7"] = dx6p3m_P7
ds_list["dx6.7m_P11"] = common.get_moni_data_list(f"{DATA_DIR}/Dx6.7m_P11/", 22, 6.25e-2, 13, 19*3600, dx6p3m_P7.sel(sec=19*3600))

ds_list["dx3.1m_P7"] = common.get_moni_data_list(f"{DATA_DIR}/Dx3.1m_P7/", 88, 3.90625e-2, 13, 19*3600, dx6p3m_P7.sel(sec=19*3600))

out_figdir="./fig_energy_budget/"
out_fext = "svg"

#---
color_list={
    "dx25m_P3": "black", "dx25m_P3_mf": "cyan", "dx25m_P7": "black", "dx25m_P7_3": "red", "dx25m_P7_0": "blue", 
    "dx25m_P7_new": "green", "dx27m_P11": "black", 
    "dx12.5m_P3": "red", "dx12.5m_P7": "red", "dx13m_P11": "red", 
    "dx6.3m_P3": "green", "dx6.3m_P7": "green", "dx6.7m_P11": "green", 
    "dx3.1m_P7": "blue",
}

linestyle_list={
    "dx25m_P3": ":", "dx25m_P7": "--",  "dx25m_P7_3": "--", "dx25m_P7_0": "--", "dx27m_P11": "-", 
    "dx25m_P7_new": "-", "dx12.5m_P3": ":", "dx12.5m_P7": "--", "dx13m_P11": "-", 
    "dx6.3m_P3": ":", "dx6.3m_P7": "--", "dx6.7m_P11": "-", 
    "dx3.1m_P7": "--",
}

exp_label_list = {
  "dx25m_P3":  "Δ=25m,P3",  
  "dx25m_P7_0":  "Δ=25m,P7,0",   
  "dx25m_P7":  "Δ=25m,P7",   
  "dx25m_P7_3":  "Δ=25m,P7,3",   
  "dx25m_P7_new":  "Δ=25m,P7,new",
  "dx27m_P11":  "Δ=27m,P11",  
  "dx12.5m_P3":  "Δ=13m,P3",
  "dx12.5m_P7":  "Δ=13m,P7",     
  'dx13m_P11': 'Δ=13m,P11', 
  'dx6.3m_P3': 'Δ=6.3m,P3',
  'dx6.3m_P7': 'Δ=6.3m,P7',  
  'dx6.7m_P11': 'Δ=6.7m,P11', 
  'dx3.1m_P7': 'Δ=3.1m,P7', 
}

if not os.path.exists(out_figdir):
    os.makedirs(out_figdir)

common.mkgraph( ds_list, "ENGK", [0.0, 1.1e11], color_list, linestyle_list, exp_label_list, 
                f"{out_figdir}/energy_budget_ENGK.{out_fext}")
for key, ds in ds_list.items():
    kmean = ds.ENGK.mean()
    print(f"{key}: {kmean.values/1e10:12.2f} x 10^10")