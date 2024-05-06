import sys
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import numpy as np

def get_numerror_data(in_numerror_logfile):
  names=['dummy1', 'tsec' ]
  for v in ["DDENS", "U", "V", "W", "DRHOT"]:
    names.append(f'L1_error_{v}')
    names.append(f'L2_error_{v}')
    names.append(f'Linf_error_{v}')
    names.append(f'Ediff_{v}')
    names.append(f'Edisp_{v}')    

  df = pd.read_csv(in_numerror_logfile, names=names, skiprows=range(0, 1), delim_whitespace=True)
  df = df.drop(columns=['dummy1']).set_index('tsec')
  ds = xr.Dataset.from_dataframe(df)
  return ds

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

def set_ax(ax, ylabel, ylim):
  ax.set_xscale("log")
  ax.set_yscale("log")
  ax.set_xlim(2.5e1, 3.5e2)
  ax.set_ylim(ylim[0], ylim[1])
  ax.tick_params(which="both", labelsize=13)
  ax.set_xlabel("$N_{e,h}\,(p+1)$", fontsize=15)
  ax.set_ylabel(ylabel, fontsize=15)
  ax.grid(which="major", color="gray", linestyle="-")
  ax.grid(which="minor", color="gray", linestyle="--")

def get_color_and_linestyle_and_label(key):
  if "P1" in key:
    color = "cyan"; lbl="p=1"
  if "P3" in key:
    color = "blue"; lbl="p=3"
  if "P7" in key:
    color = "red"; lbl="p=7"
  if "P11" in key:
    color = "green"; lbl="p=11"

  if "_dtx0.5" in key:
    lstyle = "--"; lbl=f"{lbl},C_{{rh,cs}}/2"; marker="s"
  elif "_dtx0.25" in key:
    lstyle = "-."; lbl=f"{lbl},C_{{rh,cs}}/4"; marker="D"
  else:
    lstyle = "-"; marker="o"
       
  return color, lstyle, lbl, marker

def mkgraph_num_convergence(data_list, slope_dof_list, slope_list, ylabel_list, ylim_list, outpng_path):
  print(f"mkgraph: {outpng_path} ..")
  
  fig = plt.figure(figsize=(14,6)) 
  
  ax_list = {}
  data_num = len(data_list)
  for i in range(0,data_num):
    ax_list[i] = fig.add_subplot(1, data_num, i+1)
    set_ax(ax_list[i], ylabel_list[i], ylim_list[i])

  for i in range(0,data_num):  
    for key, data in data_list[i].items():
      numerror = data
      color, lstyle, lbl, marker = get_color_and_linestyle_and_label(key)
      ax_list[i].plot( data.DOF, data, 
        label=f"${lbl}$", 
        color=color, 
        linestyle=lstyle, 
        marker=marker,         
        markeredgecolor="k", markersize=8 )
    
    for key, slope in slope_list[i].items():
      ax_list[i].plot( slope_dof_list[i], slope, 
  #      label=f"{key}", 
        color="black", 
        linestyle="--" )
  
    if i==data_num-1:
      ax_list[i].legend(bbox_to_anchor=(1.07, 1), loc='upper left', borderaxespad=0, fontsize=14)
  
  plt.subplots_adjust(wspace=0.35, hspace=0.6)
  plt.savefig(outpng_path, bbox_inches='tight')
  
def get_numerror_data_list(dir, exp_list, analysis_data):
  list_numerror = {}
  for subdir in exp_list:
    list_numerror[subdir] = get_numerror_data(f"{dir}/{subdir}/{analysis_data}")
  return list_numerror
