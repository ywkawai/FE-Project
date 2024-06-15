import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import os


def get_var1D_list(exp_name, run_s, run_e, analysis_dir, ncname):
  ds_list = []
  for runno in range(run_s, run_e+1):
      print(f'{exp_name}/run{runno}/{analysis_dir}/{ncname}.nc')
      ds = xr.open_mfdataset(f'{exp_name}/run{runno}/{analysis_dir}/{ncname}.pe*.nc', decode_times=False, combine='by_coords')
      time = ds.coords["time"].values
      print(f"time len={len(time)}")
      ds_list.append( ds.isel(time=slice(0,len(time)-1)) )
  return xr.concat(ds_list, dim="time")
      
def get_pbl_analysis(topdir, exp_list):
  pt_list = {}

  eddy_momz_flux = {}
  sgs_momz_flux = {}

  eddy_heat_flux = {}
  sgs_heat_flux = {}
  mean_heat_flux = {}
  
  for exp_dir, exp_info in exp_list.items():
    v1D_list_ = get_var1D_list(f"{topdir}/{exp_dir}", exp_info["run_sm"], exp_info["run_e"], "analysis", "v1D")
    v1D_tb_list_ = get_var1D_list(f"{topdir}/{exp_dir}", exp_info["run_sm"], exp_info["run_e"], "analysis", "tb_V1D")

    start_time = exp_info["start_time"]    
    run_sm = exp_info["run_sm"]
    run_s = exp_info["run_s"]
    run_e = exp_info["run_e"]
    integ_time_per_run = exp_info["integ_time_per_run"]
    
    time_s = start_time + (run_sm-run_s) * integ_time_per_run
    time_e = start_time + (run_e-run_s+1) * integ_time_per_run
    
    v1D_list = v1D_list_.sel(time=slice(time_s,time_e))
    v1D_tb_list = v1D_tb_list_.sel(time=slice(time_s,time_e))
    
    print(f"time:{v1D_list_.time.values}")
    print(f"time_s:{time_s}, time_e:{time_e}")
    print(v1D_list["PT"].values)
    pt_list[exp_dir] = v1D_list["PT"].mean(["time"])

    eddy_momz_flux[exp_dir] = v1D_list["MOMZ_EDDYFLX"].mean(["time"])
    sgs_momz_flux[exp_dir] = -v1D_tb_list["SGS_MOMZ_EDDYFLX"].mean(["time"])

    mean_heat_flux[exp_dir] = v1D_list["HEAT_MEANFLX"].mean(["time"])
    eddy_heat_flux[exp_dir] = v1D_list["HEAT_EDDYFLX"].mean(["time"])
    sgs_heat_flux[exp_dir] = -v1D_tb_list["SGS_HEAT_EDDYFLX"].mean(["time"])
      
  return pt_list, eddy_momz_flux, sgs_momz_flux, mean_heat_flux, eddy_heat_flux, sgs_heat_flux


def set_tick(ax, xtick_major, xtick_minor, ytick_major, ytick_minor):
  ax.set_xticks(xtick_major)
  ax.set_xticks(xtick_minor, minor=True)
  ax.xaxis.set_ticks_position('both')
  ax.set_yticks(ytick_major)
  ax.set_yticks(ytick_minor, minor=True)
  ax.yaxis.set_ticks_position('both')
  
def create_fig_vertical_structure_fillmaxmin( var_list_scaledg, vmin, vmax, tick_major, tick_minor, fig_title, exp_color_list, exp_ltype_style, exp_ltype_width, exp_label_list, out_pngname ):
  fig, ax = plt.subplots(figsize=(10.0, 8.0))
    
  for exp_name, var in var_list_scaledg.items():
    print(exp_name)
    z = var.coords["x"]
    ax.plot( var, z, 
          linestyle=exp_ltype_style[exp_name], color=exp_color_list[exp_name], label=exp_label_list[exp_name], linewidth=exp_ltype_width[exp_name])

  ax.legend(fontsize=18)
  ax.set_xlim(vmin, vmax)
  ax.set_ylim(-10, 1610)
  ax.set_ylabel("height [m]", fontsize=18) 
  ax.tick_params(which="both", labelsize=18)    
  set_tick(ax, tick_major, tick_minor, 
               np.linspace(0,1500,4), np.linspace(0,1500,16) )
  ax.set_title(fig_title, fontsize=18)
  ax.legend()
  
  plt.savefig(out_pngname)