import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

SCALERM_TMPDAT_DIR="./KT2021_SCALERM/"
CNTRL_expname = "RK8CD8_ND8Gam2e-4"
expname_list_SCALERM = [
  "RK3CD2_ND2Gam2e-4", 
  "RK4UD3", 
  "RK4CD4_ND4Gam2e-4", 
  "RK4CD4_ND8Gam2e-4", 
  "RK6UD5", 
  "RK6CD6_ND6Gam2e-4", 
  "RK8UD7", 
  "RK8CD8_ND8Gam2e-4",
]

expname_list_SCALEDG = [
  "E480P1",     
  "E480P1_MFoff",  
  "E240P3",  
  "E192P4",   
  "E160P5",  
  "E120P7", 
  "E80P11",  
]
exp_color_list = {
  "E480P1": "cyan",        
  "E480P1_MFoff": "lightskyblue",      
  "E240P3": "blue", 
  "E192P4": "red",  
  "E160P5": "green",   
  "E120P7": "goldenrod", 
  "E80P11": "black"  
}
exp_lbl_list = {
  "E480P1": "P1",      
  "E480P1_MFoff": "P1_MFoff",    
  "E240P3": "P3", 
  "E192P4": "P4",   
  "E160P5": "P5",  
  "E120P7": "P7",
  "E80P11": "P11"  
}
exp_ltype_width = {
  "E80P11": 2,  
  "E120P7": 4,
  "E160P5": 2,
  "E192P4": 4,    
  "E240P3": 2,
  "E480P1": 1,
  "E480P1_MFoff": 2,
  "reference": 1, 
}

#OUT_DIR="./compari"
OUT_DIR="./compari_FR"

def get_var_scalerm(exp_name, varname, fsuffix, varname2=None):
  vname = varname
  if varname2 != None:
    vname = varname2
    
  print(f'{SCALERM_TMPDAT_DIR}/tmp3_{exp_name}/{varname}{fsuffix}.nc')
  return xr.open_mfdataset(f'{SCALERM_TMPDAT_DIR}/tmp3_{exp_name}/{varname}{fsuffix}.nc', decode_times=False, combine='by_coords')[vname]

def get_var_scaledg(exp_name, varname, fsuffix):
#  return xr.open_mfdataset(f'../run8_{exp_name}/analysis/history{fsuffix}.pe*.nc', decode_times=False, combine='by_coords')[varname]  
  return xr.open_mfdataset(f'../run8_{exp_name}/analysis_FR/history{fsuffix}.pe*.nc', decode_times=False, combine='by_coords')[varname]


def set_tick(ax, xtick_major, xtick_minor, ytick_major, ytick_minor):
  ax.set_xticks(xtick_major)
  ax.set_xticks(xtick_minor, minor=True)
  ax.xaxis.set_ticks_position('both')
  ax.set_yticks(ytick_major)
  ax.set_yticks(ytick_minor, minor=True)
  ax.yaxis.set_ticks_position('both')

def create_fig_vertical_structure_fillmaxmin( var_list_scaledg, var_list_scalerm, vmin, vmax, tick_major, tick_minor, fig_title, exp_color_list, exp_label_list, out_pngname ):
  fig, ax = plt.subplots(figsize=(12.0, 8.0))
  
  expname0 = list(var_list_scalerm.keys())[0]
  v = var_list_scalerm[expname0]
  v_max = 1.0*v; v_min = 1.0*v
  for exp_name, var in var_list_scalerm.items():
    v_max = np.maximum(v_max, var)
    v_min = np.minimum(v_min, var)
  ax.fill_betweenx( v.coords[v.dims[0]], v_min, v_max, color="grey", alpha=0.5 )
  
  for exp_name, var in var_list_scaledg.items():
    print(exp_name)
    z = var.coords[var.dims[0]]
    ax.plot( var, z, 
          linestyle='-', color=exp_color_list[exp_name], label=exp_label_list[exp_name], linewidth=exp_ltype_width[exp_name])

  for exp_name in ["RK8CD8_ND8Gam2e-4"]:
    var = var_list_scalerm[exp_name]
    z = var.coords[var.dims[0]]
    ax.plot( var, z,   
          linestyle='--', color="black", linewidth=3, alpha=0.5) #, label="CD8ND8")

  ax.legend(fontsize=18)
  ax.set_xlim(vmin, vmax)
  ax.set_ylim(0, 1600)
#  ax.set_ylim(0, 200)  
  ax.set_ylabel("height [m]", fontsize=18) 
  ax.tick_params(which="both", labelsize=18)    
  set_tick(ax, tick_major, tick_minor, 
               np.linspace(0,1500,4), np.linspace(0,1600,17) )
  ax.set_title(fig_title, fontsize=18)
  
  plt.savefig(out_pngname)
  
  
def create_heat_flux_fig():
  tot_htflux_list_scalerm = {}
  for expname in expname_list_SCALERM:
    pt_w_prim = get_var_scalerm(expname, "PT_W_PRIM", "_XYTmean")
    sgs_htflux = get_var_scalerm(expname, "SGS_ZFLX_RHOT", "_XYTmean")
    tot_htflux_list_scalerm[expname] = pt_w_prim + sgs_htflux

  tot_htflux_list_scaledg = {}
  for expname in expname_list_SCALEDG:
    pt_w_tot = get_var_scaledg(expname, "HEAT_TOTFLX", "_v1D")
    pt_w_mean = get_var_scaledg(expname, "HEAT_MEANFLX", "_v1D")
    # pt_w_prim = get_var_scaledg(expname, "HEAT_EDDYFLX", "_v1D")
    pt_w_prim = pt_w_tot - pt_w_mean
    sgs_htflux = - get_var_scaledg(expname, "SGS_HEAT_EDDYFLX", "_tb_V1D")
    tot_htflux_list_scaledg[expname] = (pt_w_prim + sgs_htflux).mean(["time"])

  # create_fig_vertical_structure_fillmaxmin( 
  #   tot_htflux_list_scaledg, tot_htflux_list_scalerm, -30, 210, np.linspace(-0,200,5), np.linspace(-30,210,25), 
  #   "vertical heat flux [W m$^2$]", exp_color_list, exp_lbl_list, f"{OUT_DIR}/vertical_heat_flux.png" )
  create_fig_vertical_structure_fillmaxmin( 
    tot_htflux_list_scaledg, tot_htflux_list_scalerm, -30, 210, np.linspace(-0,200,5), np.linspace(-30,210,25), 
    "vertical heat flux [W m$^2$]", exp_color_list, exp_lbl_list, f"{OUT_DIR}/vertical_heat_flux_zoom.png" )

def create_w_prim2_fig():
  w_w_prim_list_scalerm = {}
  for expname in expname_list_SCALERM:
    print(expname)
    w_w_prim_list_scalerm[expname] = get_var_scalerm(expname, "W_PRIM2", "_XYTmean").mean(["time"])

  w_w_prim_list_scaledg = {}
  for expname in expname_list_SCALEDG:
    w_w_tot = get_var_scaledg(expname, "MOMZ_TOTFLX", "_v1D")
    w_w_mean = get_var_scaledg(expname, "MOMZ_MEANFLX", "_v1D")
    dens_mean = get_var_scaledg(expname, "DENS", "_v1D")    
    w_w_prim = ( w_w_tot - w_w_mean ) / dens_mean
    w_w_prim_list_scaledg[expname] = w_w_prim.mean(["time"])

  create_fig_vertical_structure_fillmaxmin( 
    w_w_prim_list_scaledg, w_w_prim_list_scalerm, 0.0,2.0, np.linspace(0.0,2.0,5), np.linspace(0.0,2.0,21), 
    "W_PRIM2 [m$^2$/s$^2$]", exp_color_list, exp_lbl_list, f"{OUT_DIR}/w_prim2.png" )

def create_w_skewness_fig():
  w_skew_list_scalerm = {}
  for expname in expname_list_SCALERM:    
    w_prim2 = get_var_scalerm(expname, "W_PRIM2", "_XYTmean")
    w_prim3 = get_var_scalerm(expname, "W_PRIM3", "_XYTmean")
    w_skew_list_scalerm[expname] =  (w_prim3/w_prim2**1.5).mean(["time"])

  w_skew_list_scaledg = {}
  for expname in expname_list_SCALEDG:
    w_prim3 = get_var_scaledg(expname, "W_PRIM3", "_v1D") 
    w_w_tot = get_var_scaledg(expname, "MOMZ_TOTFLX", "_v1D")
    w_w_mean = get_var_scaledg(expname, "MOMZ_MEANFLX", "_v1D")  
    dens_mean = get_var_scaledg(expname, "DENS", "_v1D")        
    w_prim2 = ( w_w_tot - w_w_mean ) / dens_mean
    w_skew_list_scaledg[expname] =  (w_prim3/w_prim2**1.5).mean(["time"])

  create_fig_vertical_structure_fillmaxmin( 
    w_skew_list_scaledg, w_skew_list_scalerm, -0.5,1.5, np.linspace(-0.5,1.5,5), np.linspace(-0.5,1.5,21), 
    "Skewness", exp_color_list, exp_lbl_list, f"{OUT_DIR}/w_skewness.png" )

def create_pt_fig():
  pt_list_scalerm = {}
  for expname in expname_list_SCALERM:
    pt_list_scalerm[expname] = get_var_scalerm(expname, "RHOTovDENS", "_XYTmean", "PT")

  pt_list_scaledg = {}
  for expname in expname_list_SCALEDG:
    pt_list_scaledg[expname] = get_var_scaledg(expname, "PT", "_v1D").mean(["time"])

  create_fig_vertical_structure_fillmaxmin( 
    pt_list_scaledg, pt_list_scalerm, 303.5,305.5, np.linspace(303.5,305.5,5), np.linspace(303.5,305.5,21), 
    "potential temperature [K]", exp_color_list, exp_lbl_list, f"{OUT_DIR}/PotTemp.png" )

#-----  
#create_pt_fig()
#create_heat_flux_fig()
#create_w_prim2_fig()
create_w_skewness_fig()