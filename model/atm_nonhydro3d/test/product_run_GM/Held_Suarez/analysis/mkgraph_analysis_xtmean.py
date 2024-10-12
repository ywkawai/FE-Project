import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import matplotlib.animation as animation
import os

EXP_TOP_DIR = "./rhot_hevi"
EXP_DIR_LIST = {
#  "Eh12Ez8P3", 
#  "Eh24Ez16P3",     
#    "Eh48Ez32P3"
   "Eh6Ez4P7",  
  # "Eh12Ez8P7", 
  # "Eh24Ez16P7", 
  # "Eh48Ez32P7",
#  "Eh4Ez3P11", "Eh8Ez6P11", 
#  "Eh16Ez12P11"
}

ANALYSIS_OUT_DIR="analysis_out"

#--------------------------------------
def v_1Daxis_fmt(tick_val, pos):
  val = int(tick_val/100)
  return f'{val}'

def set_fig_Yp_axis(ax):
  ax.tick_params(labelsize=11, length=2)
  ax.set_xlabel('latitude [degrees]', fontsize=14) 
  ax.set_xlim([-89,89]) 
  ax.xaxis.set_major_locator(tick.MultipleLocator(30.0))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(10.0))
  ax.set_ylabel('pressure [hPa]', fontsize=14)
  ax.invert_yaxis()
  ax.yaxis.set_major_locator(tick.MultipleLocator(100e2))
  ax.yaxis.set_major_formatter(tick.FuncFormatter(v_1Daxis_fmt))

def plot_var_yz(v, fig_title, vmin, vmax, vint, interval, cbar_int, pdf_name):
  print(f"plot: {pdf_name}")

  y = v.coords["lat"]
  z = v.coords["p"]

  Y, Z = np.meshgrid(y,z)
  fig = plt.figure(figsize=(9,4.2)) 
  ax = fig.add_subplot(1,1,1)
  set_fig_Yp_axis(ax)
  ax.set_title(fig_title, fontsize=14)

  lv = np.linspace(vmin,vmax,int((vmax-vmin)/vint)+1)      
  pcm = ax.contourf(Y, Z, v, levels=lv, cmap='jet', extend='both')
  cont = ax.contour(Y, Z, v, levels=np.arange(vmin, vmax+interval, interval), colors=['black'])
  cont.clabel(fmt='%2.0f', fontsize=12)
  
  fmt = tick.ScalarFormatter(useMathText=True)
  #fmt.set_powerlimits((0,0))
  cbar = plt.colorbar(pcm, aspect=40.0, extend='both', orientation='vertical', shrink=0.8, format=fmt, pad=0.02)
  cbar.ax.tick_params(labelsize=11)
  cbar.ax.yaxis.get_offset_text().set_fontsize(11)
  cbar.ax.yaxis.set_tick_params(pad=2)     
  cbar.set_ticks(np.arange(vmin, vmax+1e-10, cbar_int))    

  plt.savefig(pdf_name)

def mean_lon( var ):
  lon = var["lon"]
  lon_weight = [ 0.166666666666667, 0.833333333333333, 0.833333333333333, 0.166666666666667 ]
  nelem_lon = int(len(lon)/4)

  df = pd.DataFrame({ 'lon': lon.values,  'lon_intweight': lon_weight  *  nelem_lon } ).set_index("lon")
  ds = df.to_xarray()
  var_mean_lon = ( var * ds["lon_intweight"] ).sum("lon") / (2.0 * nelem_lon )
  var_mean_lon.name = var.name
  return var_mean_lon

def analysis_xtmean(exp_dir, out_dir):
  vars_list = ["Umet", "Vmet", "T", "merid_eddy_momflx", "merid_eddy_hflx", "eddy_kinetic_energy", "eddy_temp_variance"]
  vars = {}
  for varname in vars_list:
    vars[varname] = xr.open_mfdataset(f"{out_dir}/tmp_data/{varname}.nc", decode_times=True, use_cftime=True, combine='by_coords')[varname]    

  out_fig_dir = f"{out_dir}/figs"
  os.makedirs(out_fig_dir, exist_ok=True)      
  plot_var_yz(vars["Umet"].mean("time"), 'Zonal wind [m/s]', -36.0, 36.0, 1.0, 4.0, 8, f"{out_fig_dir}/Umet_xtmean.pdf")
  plot_var_yz(vars["Vmet"].mean("time"), 'Meridonal wind [m/s]', -4.5, 4.5, 0.25, 0.25,  0.5, f"{out_fig_dir}/Vmet_xtmean.pdf")
  plot_var_yz(vars["T"].mean("time"), 'Temperature [m/s]', 180.0, 310.0, 1.25, 5.0,  10.0, f"{out_fig_dir}/T_xtmean.pdf")
  plot_var_yz(vars["merid_eddy_momflx"].mean("time"), 'Merid. Eddy Flux of Zonal Momentum [m2/s2]', -70.0, 70.0, 2.5, 10.0, 10.0, f"{out_fig_dir}/Merid_Eddy_momflx_xtmean.pdf")
  plot_var_yz(vars["merid_eddy_hflx"].mean("time"), 'Merid. Eddy Flux of Temperature [K m/s]', -25.0, 25.0, 0.625, 2.5,  5.0, f"{out_fig_dir}/Merid_Eddy_hflx_xtmean.pdf")
  plot_var_yz(vars["eddy_kinetic_energy"].mean("time"), 'Eddy Kinetic Energy[m2/s2]', 0.0, 420.0, 6.0, 30.0, 60,  f"{out_fig_dir}/Eddy_Kinetic_Energy_xtmean.pdf")
  plot_var_yz(vars["eddy_temp_variance"].mean("time"), 'Eddy Temperature Variance[K2]', 0.0, 45.0, 1.0, 5.0,  5.0, f"{out_fig_dir}/Eddy_Temp_variance_xtmean.pdf")

for exp_dir in EXP_DIR_LIST:
  analysis_xtmean(f"{EXP_TOP_DIR}/{exp_dir}", f"{ANALYSIS_OUT_DIR}/{exp_dir}")