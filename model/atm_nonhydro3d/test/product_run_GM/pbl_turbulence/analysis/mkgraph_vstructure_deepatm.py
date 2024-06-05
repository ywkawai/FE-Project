import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
import os

expname_list_SCALEDG_noSA = [
  "Eh128Ez64P3_deepatm",     
  "Eh100Ez52P4_deepatm",       
  "Eh64Ez34P7_deepatm", 
]

expname_list_SCALEDG = [
  "Eh128Ez64P3",     
  "Eh100Ez52P4",       
  "Eh64Ez34P7", 
]
exp_color_list = {
  "Eh100Ez52P4_deepatm": "red",  
  "Eh128Ez64P3_deepatm": "blue",     
#  "Eh64Ez34P7_deepatm": "goldenrod",   
  "Eh64Ez34P7_deepatm": "goldenrod", 
}
exp_lbl_list = {
  "Eh128Ez64P3_deepatm": "$p=3$",  
  "Eh100Ez52P4_deepatm": "$p=4$",      
  "Eh64Ez34P7_deepatm": "$p=7$",
}
exp_ltype_width = {
  "Eh64Ez34P7_deepatm": 4,    
  "Eh100Ez52P4_deepatm": 4,    
  "Eh128Ez64P3_deepatm": 2,
}

#OUT_DIR="./compari"
OUT_DIR="./analysis_out/vprofile_noSA"
EXP_TOPDIR="rp3.4km"

#---------------------
def get_var_scaledg_KT2023(ncpath, varname):
  print(ncpath)
  return xr.open_mfdataset(ncpath, decode_times=False, combine='by_coords')[varname]

def get_var_scaledg(exp_name, varname, fsuffix):
  ncpath = f'{EXP_TOPDIR}/{exp_name}/run8/analysis/history{fsuffix}.pe*.nc'
  print(ncpath)
  return xr.open_mfdataset(ncpath, decode_times=False, combine='by_coords')[varname]

def set_tick(ax, xtick_major, xtick_minor, ytick_major, ytick_minor):
  ax.set_xticks(xtick_major)
  ax.set_xticks(xtick_minor, minor=True)
  ax.xaxis.set_ticks_position('both')
  ax.set_yticks(ytick_major)
  ax.set_yticks(ytick_minor, minor=True)
  ax.yaxis.set_ticks_position('both')
  ax.tick_params(axis="both", which="major", length=7)
  ax.tick_params(axis="both", which="minor", length=5)


def remove_duplicate_pos(v):
  x = v.coords["x"].values
  x_uniq = pd.unique(x)
  x_ind = []
  for xp in x_uniq:
    x_ind.append(np.where(x==xp)[0][0])
  return v.isel(x=x_ind)

def create_fig_vertical_structure_fillmaxmin( var_list_scaledg, var_list_scaledg_KT2023, vmin, vmax, tick_major, tick_minor, fig_title, exp_color_list, exp_label_list, out_pngname ):
  fig, ax = plt.subplots(figsize=(12.0, 8.0))
  
  expname0 = list(var_list_scaledg_KT2023.keys())[0]
  v_ = var_list_scaledg_KT2023[expname0]
  v = remove_duplicate_pos(v_)
#  print(v.coords["x"])  
  v_max = 1.0*v; v_min = 1.0*v
  
  for exp_name, var_ in var_list_scaledg_KT2023.items():
    var = remove_duplicate_pos(var_).interp(x=v.coords["x"].values)
    print(var.coords["x"])
    v_max = np.maximum(v_max, var)
    v_min = np.minimum(v_min, var)
#  print(v_max)
  ax.fill_betweenx( v.coords[v.dims[0]], v_min, v_max, color="gray", alpha=0.5 )
  
  for exp_name, var in var_list_scaledg.items():
#    print(exp_name)
    z = var.coords[var.dims[0]]
    ax.plot( var, z, 
          linestyle='-', color=exp_color_list[exp_name], label=exp_label_list[exp_name], linewidth=exp_ltype_width[exp_name])

  ax.legend(fontsize=24)
  ax.set_xlim(vmin, vmax)
  ax.set_ylim(0, 1600)
#  ax.set_ylim(0, 200)  
  ax.set_ylabel("height [m]", fontsize=24) 
  ax.tick_params(which="both", labelsize=22)    
  set_tick(ax, tick_major, tick_minor, 
               np.linspace(0,1500,4), np.linspace(0,1600,17) )
  ax.set_title(fig_title, fontsize=20)
  
  plt.savefig(out_pngname)
  
  
#-
def cal_w_prim2_dg(expname_list):
  w_w_prim_list_scaledg = {}
  for expname in expname_list:
    w_w_tot = get_var_scaledg(expname, "MOMZ_TOTFLX", "_v1D")
    w_w_mean = get_var_scaledg(expname, "MOMZ_MEANFLX", "_v1D")
    dens_mean = get_var_scaledg(expname, "DENS", "_v1D")    
    w_w_prim = ( w_w_tot - w_w_mean ) / dens_mean
    w_w_prim_list_scaledg[expname] = w_w_prim.mean(["time"])
  
  return w_w_prim_list_scaledg
  
def create_w_prim2_fig():
  w_w_prim_list_scaledg = cal_w_prim2_dg( expname_list_SCALEDG )   
  w_w_prim_list_scaledg_noSA = cal_w_prim2_dg( expname_list_SCALEDG_noSA )  
  create_fig_vertical_structure_fillmaxmin( 
    w_w_prim_list_scaledg_noSA, w_w_prim_list_scaledg, 0.0,2.0, np.linspace(0.0,2.0,5), np.linspace(0.0,2.0,21), 
    "W_PRIM2 [m$^2$/s$^2$]", exp_color_list, exp_lbl_list, f"{OUT_DIR}/w_prim2.pdf" )

#-
def cal_w_skewness_dg(expname_list):
  w_skew_list_scaledg = {}
  for expname in expname_list:
    w_prim3 = get_var_scaledg(expname, "W_PRIM3", "_v1D") 
    w_w_tot = get_var_scaledg(expname, "MOMZ_TOTFLX", "_v1D")
    w_w_mean = get_var_scaledg(expname, "MOMZ_MEANFLX", "_v1D")  
    dens_mean = get_var_scaledg(expname, "DENS", "_v1D")        
    w_prim2 = ( w_w_tot - w_w_mean ) / dens_mean
    w_skew_list_scaledg[expname] =  (w_prim3/w_prim2**1.5).mean(["time"])
  return w_skew_list_scaledg
  
def create_w_skewness_fig():
  w_skew_list_scaledg = cal_w_skewness_dg( expname_list_SCALEDG ) 
  w_skew_list_scaledg_noSA = cal_w_skewness_dg( expname_list_SCALEDG_noSA ) 
  create_fig_vertical_structure_fillmaxmin( 
    w_skew_list_scaledg_noSA, w_skew_list_scaledg, -0.5,1.5, np.linspace(-0.5,1.5,5), np.linspace(-0.5,1.5,21), 
    "Skewness", exp_color_list, exp_lbl_list, f"{OUT_DIR}/w_skewness.pdf" )

#-
def cal_heat_flux_dg(expname_list):
  tot_htflux_list_scaledg = {}
  for expname in expname_list:
    pt_w_tot = get_var_scaledg(expname, "HEAT_TOTFLX", "_v1D")
    pt_w_mean = get_var_scaledg(expname, "HEAT_MEANFLX", "_v1D")
    # pt_w_prim = get_var_scaledg(expname, "HEAT_EDDYFLX", "_v1D")
    pt_w_prim = pt_w_tot - pt_w_mean
    sgs_htflux = - get_var_scaledg(expname, "SGS_HEAT_EDDYFLX", "_tb_V1D")
    tot_htflux_list_scaledg[expname] = (pt_w_prim + sgs_htflux).mean(["time"])
  return tot_htflux_list_scaledg
  
def create_heat_flux_fig():
  tot_htflux_list_scaledg = cal_heat_flux_dg(expname_list_SCALEDG)
  tot_htflux_list_scaledg_noSA = cal_heat_flux_dg(expname_list_SCALEDG_noSA)
  
  create_fig_vertical_structure_fillmaxmin( 
    tot_htflux_list_scaledg_noSA, tot_htflux_list_scaledg, -30, 210, np.linspace(-0,200,5), np.linspace(-30,210,25), 
    "vertical heat flux [W m$^2$]", exp_color_list, exp_lbl_list, f"{OUT_DIR}/vertical_heat_flux_zoom.pdf" )

#-
def cal_pt_dg(expname_list):
  pt_list_scaledg = {}
  for expname in expname_list:
    pt_list_scaledg[expname] = get_var_scaledg(expname, "PT", "_v1D").mean(["time"])
  return pt_list_scaledg

def create_pt_fig():
  pt_list_scaledg = cal_pt_dg(expname_list_SCALEDG)
  pt_list_scaledg_noSA = cal_pt_dg(expname_list_SCALEDG_noSA)
  create_fig_vertical_structure_fillmaxmin( 
    pt_list_scaledg_noSA, pt_list_scaledg, 303.0,305.5, np.arange(303.0,305.5,0.5), np.arange(303.0,305.5,0.1), 
    "potential temperature [K]", exp_color_list, exp_lbl_list, f"{OUT_DIR}/PotTemp.pdf" )
  
#-----
os.makedirs(OUT_DIR, exist_ok=True)
create_pt_fig()
create_heat_flux_fig()
create_w_prim2_fig()
create_w_skewness_fig()