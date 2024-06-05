import sys
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import numpy as np
import os

OUT_DIR="analysis_out/ps_l2error"
VAR_LIST = ["SFCPRES"]


# p=3
P3_exp_list = ["Eh10Ez8P3", "Eh20Ez8P3", "Eh40Ez8P3", "Eh80Ez8P3"]
P3_exp_label = {
  "Eh10Ez8P3": "$\Delta_{h,eq}\sim 250$ km", 
  "Eh20Ez8P3": "$\Delta_{h,eq}\sim 125$ km", 
  "Eh40Ez8P3": "$\Delta_{h,eq}\sim 63$ km", 
  "Eh80Ez8P3": "$\Delta_{h,eq}\sim 32$ km",   
}
P3_exp_lstyle = {
  "Eh10Ez8P3": "-", 
  "Eh20Ez8P3": "--", 
  "Eh40Ez8P3": ":", 
  "Eh80Ez8P3": "-.",   
}


# p=7
P7_exp_list = ["Eh5Ez4P7", "Eh10Ez4P7", "Eh20Ez4P7", "Eh40Ez4P7"]
P7_exp_label = {
  "Eh5Ez4P7": "$\Delta_{h,eq}\sim 250$ km", 
  "Eh10Ez4P7": "$\Delta_{h,eq}\sim 125$ km", 
  "Eh20Ez4P7": "$\Delta_{h,eq}\sim 63$ km", 
  "Eh40Ez4P7": "$\Delta_{h,eq}\sim 32$ km",   
}
P7_exp_lstyle = {
  "Eh5Ez4P7": "-", 
  "Eh10Ez4P7": "--", 
  "Eh20Ez4P7": ":", 
  "Eh40Ez4P7": "-.",   
}

# p=11
P11_exp_list = ["Eh3Ez3P11", "Eh6Ez3P11", "Eh12Ez3P11", "Eh24Ez3P11"]
P11_exp_label = {
  "Eh3Ez3P11": "$\Delta_{h,eq}\sim 278$ km", 
  "Eh6Ez3P11": "$\Delta_{h,eq}\sim 139$ km", 
  "Eh12Ez3P11": "$\Delta_{h,eq}\sim 69$ km", 
  "Eh24Ez3P11": "$\Delta_{h,eq}\sim 35$ km",   
}
P11_exp_lstyle = {
  "Eh3Ez3P11": "-", 
  "Eh6Ez3P11": "--", 
  "Eh12Ez3P11": ":", 
  "Eh24Ez3P11": "-.",   
}

#-----------------------------
def get_numerror_data(in_numerror_logfile):
  names=['dummy1', 'tsec' ]
  for v in VAR_LIST:
    names.append(f'L1_error_{v}')
    names.append(f'L2_error_{v}')
    names.append(f'Linf_error_{v}')
    names.append(f'Ediff_{v}')
    names.append(f'Edisp_{v}')    

  df = pd.read_csv(in_numerror_logfile, names=names, skiprows=range(0, 1), delim_whitespace=True)
  df.loc[:,"day"] = df.loc[:,"tsec"] / 86400.0   
  df = df.drop(columns=['dummy1']).set_index('day')
  ds = xr.Dataset.from_dataframe(df)
  return ds

def get_numerror_JM06_urange(in_csvfile):
  df = pd.read_csv(in_csvfile, names=["day", "L2_error_SFCPRES"])
  df = df.set_index('day')  
  return xr.Dataset.from_dataframe(df)  

def resoldep_numerrordata_sel_time(varname, time, numerror_data_list, dof_list):
  
  dof = np.array(list(dof_list.values()))
  data_size = len(dof)

  l1_error = xr.DataArray(np.zeros(data_size), dims=['DOF'], 
                  coords={'DOF':dof}, name=f'L1_error_{varname}')
  l2_error = xr.DataArray(np.zeros(data_size), dims=['DOF'], 
                  coords={'DOF':dof}, name=f'L2_error_{varname}')
  linf_error = xr.DataArray(np.zeros(data_size), dims=['DOF'], 
                  coords={'DOF':dof}, name=f'Linf_error_{varname}')  
  ediff = xr.DataArray(np.zeros(data_size), dims=['DOF'], 
                  coords={'DOF':dof}, name=f'Ediff_{varname}')  
  edisp = xr.DataArray(np.zeros(data_size), dims=['DOF'], 
                  coords={'DOF':dof}, name=f'Edisp_{varname}')  
  
  i = 0
  for key, numerror_data in numerror_data_list.items():
    l1_error[i] = numerror_data.sel(tsec=time)[f'L1_error_{varname}'].values    
    l2_error[i] = numerror_data.sel(tsec=time)[f'L2_error_{varname}'].values  
    linf_error[i] = numerror_data.sel(tsec=time)[f'Linf_error_{varname}'].values    
    ediff[i] = numerror_data.sel(tsec=time)[f'Ediff_{varname}'].values    
    edisp[i] = np.abs( numerror_data.sel(tsec=time)[f'Edisp_{varname}'].values )
    
    i = i + 1
    
#  return l1_error, l2_error, linf_error, ediff, edisp
  return l1_error, l2_error, linf_error

def merge_ds(ds_list):
  ds_list_ = []
  i = 0
  for ds in ds_list:
    time = ds["day"].values
    tlen = len(time)  
    if i > 0:
      ds = ds.isel(day=slice(1,tlen))
    
    ds_list_.append(ds)
    i = i+1
  
  return xr.concat(ds_list_, "day")

def get_numerror_merged_data( resol_dir_list, varname ):
    l2error = {}
    for subdir in resol_dir_list:
        numerror_tmp_data = []
        for runno in [1, 2]:
            ds = get_numerror_data(f"rhot_hevi/{subdir}_{runno}/analysis/NUMERROR_SFC_LOG.peall")
            numerror_tmp_data.append( ds[varname] )
        l2error[subdir] = merge_ds(numerror_tmp_data)

    return l2error

def mkgraph(var_list, var_urange, 
            label_list, line_color, line_style_list, 
            png_out):
  fig = plt.figure(figsize=(4,7))
  ax = fig.add_subplot(111)
  ax.set_yscale("log")
  ax.set_xlim([0,1e1])  
  ax.set_ylim([1e-6,2e1])

  zero = var_urange["L2_error_SFCPRES"].values * 0.0
  ax.fill_between(var_urange.day, var_urange["L2_error_SFCPRES"].values, zero+1e-6, 
                  facecolor='lightgray', alpha=0.3)
  
  for key, var in var_list.items():
    print(var)
    ax.plot(var.day, var.values / 1e2, 
            color=line_color, 
            linestyle=line_style_list[key], 
            label=label_list[key],
            linewidth=2 )
  ax.set_ylabel("L2 error (Ps) [hPa]", fontsize=14)
  ax.set_xlabel("time [day]", fontsize=14)
  ax.tick_params(axis="both", labelsize=12)
  ax.legend(fontsize=12)
  
  plt.savefig(png_out)

#------------------------------

os.makedirs(OUT_DIR, exist_ok=True)

l2error_urange = get_numerror_JM06_urange("analysis/JM06_uncertain_range.csv")

l2error_P3 = get_numerror_merged_data( P3_exp_list, "L2_error_SFCPRES" )
mkgraph( l2error_P3, l2error_urange.sel(day=slice(0,11)),
        P3_exp_label, "blue", P3_exp_lstyle, 
        f"{OUT_DIR}/ps_l2error_P3.pdf" )

l2error_P7 = get_numerror_merged_data( P7_exp_list, "L2_error_SFCPRES" )
mkgraph( l2error_P7, l2error_urange.sel(day=slice(0,11)),
        P7_exp_label, "red", P7_exp_lstyle, 
        f"{OUT_DIR}/ps_l2error_P7.pdf" )

l2error_P11 = get_numerror_merged_data( P11_exp_list, "L2_error_SFCPRES" )
mkgraph( l2error_P11, l2error_urange.sel(day=slice(0,11)),
        P11_exp_label, "green", P11_exp_lstyle, 
        f"{OUT_DIR}/ps_l2error_P11.pdf" )
