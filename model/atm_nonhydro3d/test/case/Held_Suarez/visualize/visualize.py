import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as tick

run_dir_list = [ "run1", "run2" ]
dist_dir = "./visualize"

#------------------

def v_1Daxis_fmt(tick_val, pos):
  val = int(tick_val/100)
  return f'{val}'

def set_fig_Yp_axis(ax):
  ax.tick_params(labelsize=15, length=8)
  ax.set_xlabel('latitude [degrees]', fontsize=20)  
  ax.xaxis.set_major_locator(tick.MultipleLocator(30.0))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(10.0))
  ax.set_ylabel('pressure [hPa]', fontsize=20)
  ax.invert_yaxis()
  ax.yaxis.set_major_locator(tick.MultipleLocator(100e2))
  ax.yaxis.set_major_formatter(tick.FuncFormatter(v_1Daxis_fmt))

def plot_var_yz(v, fig_title, vmin, vmax, interval, png_name):
  y = v.coords["lat"]
  z = v.coords["p"]

  Y, Z = np.meshgrid(y,z)
  fig = plt.figure(figsize=(16,8)) 
  ax = fig.add_subplot(1,1,1)
  set_fig_Yp_axis(ax)
  ax.set_title(fig_title, fontsize=22)

  print(f"output: {png_name} ..")
  pcm = ax.pcolormesh(Y, Z, v.values, vmin=vmin, vmax=vmax, cmap='jet')
  cont = ax.contour(Y, Z, v.values, levels=np.arange(vmin, vmax+interval, interval), colors=['black'])
  cont.clabel(fmt='%2.0f', fontsize=12)

  fmt = tick.ScalarFormatter(useMathText=True)
  cbar = plt.colorbar(pcm, aspect=40.0, extend='both', shrink=0.8, format=fmt)
  cbar.ax.tick_params(labelsize=18)

  plt.savefig(png_name)


def mean_lon( var ):
  lon = var["lon"]
  lon_weight = [ 0.166666666666667, 0.833333333333333, 0.833333333333333, 0.166666666666667 ]
  nelem_lon = int(len(lon)/4)

  df = pd.DataFrame({ 'lon': lon.values,  'lon_intweight': lon_weight  *  nelem_lon } ).set_index("lon")
  ds = df.to_xarray()
  var_mean_lon = ( var * ds["lon_intweight"] ).sum("lon") / (2.0 * nelem_lon )
  var_mean_lon.name = var.name
  return var_mean_lon

  #----------------------------------------------------

vars_list = {"Umet":[], "Vmet":[], "W": [], "T":[], "PRES":[], "merid_eddy_momflx":[], 
  "PT":[], "merid_eddy_hflx":[], "eddy_kinetic_energy":[], "eddy_temp_variance":[] }
vars = {}

for n, run_dir in enumerate(run_dir_list):
  dir = f"{run_dir}/outdata_p"  
  nc = xr.open_mfdataset(f"{dir}/history.pe000*.nc", decode_times=True, combine='by_coords')
  time = nc["time"]

  if n > 0: 
    nc = nc.isel(time=slice(1,len(time.values)))   
  
  for vname in ["Umet", "Vmet", "W", "T", "PRES"]:
    vars_list[vname].append( nc[vname] )

  up = nc["Umet"] - mean_lon(nc["Umet"])
  vp = nc["Vmet"] - mean_lon(nc["Vmet"])
  tp = nc["T"] - mean_lon(nc["T"])

  up_vp = ( vp * up )
  up_vp.name = "merid_eddy_momflx"
  vars_list["merid_eddy_momflx"].append( up_vp )

  vp_tp = ( vp * tp )
  vp_tp.name = "merid_eddy_hflx"
  vars_list["merid_eddy_hflx"].append( vp_tp )

  eke = ( 0.5*(up**2 + vp**2) )
  eke.name = "eddy_kinetic_energy"
  vars_list["eddy_kinetic_energy"].append( eke )  

  tp_tp = ( tp * tp )
  tp_tp.name = "eddy_temp_variance"
  vars_list["eddy_temp_variance"].append( tp_tp )  

  pt = nc["T"] * ( 1e5 / nc["PRES"] )**(287.0/1004.0)
  pt.name = "PT"
  vars_list["PT"].append( pt )  

for key, vlist in vars_list.items():
  vars[key] = xr.concat(vlist, "time")

plot_var_yz(mean_lon(vars["Umet"]).mean(["time"]), 'Zonal wind [m/s]', -36.0, 36.0, 4.0,  f"{dist_dir}/Umet_xtmean.png")
plot_var_yz(mean_lon(vars["Vmet"]).mean(["time"]), 'Meridonal wind [m/s]', -4.25, 4.25, 0.25,  f"{dist_dir}/Vmet_xtmean.png")
plot_var_yz(mean_lon(vars["W"]).mean(["time"]), 'Vertical wind [m/s]', -4e-3, 4e-3, 4e-4,  f"{dist_dir}/W_xtmean.png")
plot_var_yz(mean_lon(vars["T"]).mean(["time"]), 'Temperature [m/s]', 185.0, 305.0, 5.0,  f"{dist_dir}/T_xtmean.png")
plot_var_yz(mean_lon(vars["merid_eddy_momflx"]).mean(["time"]), 'Merid. Eddy Flux of Zonal Momentum [m2/s2]', -80.0, 80.0, 10.0,  f"{dist_dir}/Merid_Eddy_momflx_xtmean.png")
plot_var_yz(mean_lon(vars["merid_eddy_hflx"]).mean(["time"]), 'Merid. Eddy Flux of Temperature [K m/s]', -24.0, 24.0, 3.0,  f"{dist_dir}/Merid_Eddy_hflx_xtmean.png")
plot_var_yz(mean_lon(vars["eddy_kinetic_energy"]).mean(["time"]), 'Eddy Kinetic Energy[m2/s2]', 0.0, 440.0, 40.0,  f"{dist_dir}/Eddy_Kinetic_Energy_xtmean.png")
plot_var_yz(mean_lon(vars["eddy_temp_variance"]).mean(["time"]), 'Eddy Temperature Variance[K2]', 0.0, 48.0, 4.0,  f"{dist_dir}/Eddy_Temp_variance_xtmean.png")
