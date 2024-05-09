import xarray as xr
import numpy as np
import mkgraph_common as mkgraph
import matplotlib.ticker as tick
import matplotlib.pyplot as plt
import os

target_dir="./rhot_heve/"
out_dir="./analysis_out/h25m/num_error/"
time_list = [1800, 3600, 5400, 7200]

#--------------

class analysis_numconv:
  def __init__(self, var_ref_tmp, pord=7) -> None:
    self.node_num = pord+1
    if pord==7:
        w1d = [0.0357,0.210,0.341,0.412,0.412,0.341,0.21,0.0357]
    elif pord==3:
        w1d = [0.16666,0.83333,0.83333,0.16666]
        
    w_h1d = np.array(w1d*int(var_ref_tmp.shape[1]/self.node_num))
    w_v1d = np.array(w1d*int(var_ref_tmp.shape[0]/self.node_num))
    print(var_ref_tmp.shape)
    print(w_h1d.shape)    
    print(w_v1d.shape)
    w = np.tile(w_h1d, (var_ref_tmp.shape[0],1))
    self.w_3D = w * np.transpose( np.tile(w_v1d, (var_ref_tmp.shape[1],1)) )
    self.fac=0.5*self.node_num
  def eval_numerror_conv( self, resol_list, var_list, var_ref, error_type="l1"):
      numerror_list_ = {}
      w_3D_ = self.w_3D * self.fac**2
      for resol in resol_list:
          if error_type=="l0":
              numerror = (w_3D_*(var_list[resol]-var_ref)).coarsen(x=self.node_num,z=self.node_num).mean()        
          elif error_type=="l1":
              numerror = (w_3D_*np.abs(var_list[resol]-var_ref)).coarsen(x=self.node_num,z=self.node_num).mean()
          elif error_type=="l2":
              numerror = np.sqrt((w_3D_*(var_list[resol]-var_ref)**2).coarsen(x=self.node_num,z=self.node_num).mean())
          numerror_list_[resol] = numerror
      return numerror_list_

#------------------

def get_fpathlist(
  dir1,ftype,domainlabel,timelabel,
  PRC_NUM_X,PRC_NUM_Y, 
  PRC_NUM_X_s,PRC_NUM_X_e,PRC_NUM_Y_s,PRC_NUM_Y_e):
  prc_offset = 0#(panelID-1)*PRC_NUM_X*PRC_NUM_Y  
  return [ [ dir1 + "{0}{1}{2}.pe{3:06d}.nc".format( ftype, domainlabel, timelabel, prc_offset + i + j*PRC_NUM_X ) for i in range(PRC_NUM_X_s,PRC_NUM_X_e) ] for j in range(PRC_NUM_Y_s,PRC_NUM_Y_e) ]

def merge_hist(fpath, dim=["y","x"]):
    return xr.open_mfdataset(fpath, decode_times=False,combine="nested", concat_dim=dim)


def hori_1Daxis0_fmt(tick_val, pos):
  val = int(tick_val)
  return f'{val}'
def hori_1Daxis_fmt(tick_val, pos):
  val = int(tick_val/1e3)
#  print(tick_val*Radius+120e3)  
  return f'{val}'
def vert_1Daxis_fmt(tick_val, pos):
  val = int((tick_val-6e4)/1e3)
#  print(val)
  return f'{val}'
def set_fig_XZ_axis(ax, xlabel):
  ax.tick_params(labelsize=18, length=10)  
#  print("ylabel:", xlabel)
  ax.xaxis.set_major_locator(tick.MultipleLocator(10000.0))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(2500.0))    
  ax.xaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  if (xlabel):  
    ax.set_xlabel('x [km]', fontsize=18)  
  ax.yaxis.set_major_locator(tick.MultipleLocator(2000.0))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(500.0))  
  ax.yaxis.set_major_formatter(tick.FuncFormatter(vert_1Daxis_fmt))
  ax.set_ylabel('height [km]', fontsize=18)

#------------------

fpath = get_fpathlist(target_dir,"Eh192Ez30P3/outdata_compari/history","","",64,1, 0,64,0,1)
df = merge_hist(fpath)
fpath_hi2 = get_fpathlist(target_dir,"Eh384Ez60P3/outdata_compari/history","","",64,1, 0,64,0,1)
df_hi2 = merge_hist(fpath_hi2)
fpath_hi4 = get_fpathlist(target_dir,"Eh1536Ez240P3/outdata_compari/history","","",64,1, 0,64,0,1)
df_hi4 = merge_hist(fpath_hi4)

target_dir_ref = "../mountain_wave_analytic_metric/rhot_heve/"
fpath_ref = get_fpathlist(target_dir_ref,"Eh384Ez60P7/outdata_compari/history","","",64,1, 0,64,0,1)
df_ref = merge_hist(fpath_ref)


w = df["W"].isel(y=0,time=12)
w_hi2 = df_hi2["W"].isel(y=0,time=12)
w_hi4 = df_hi4["W"].isel(y=0,time=12)
w_ref = df_ref["W"].isel(y=0,time=12)
numconv = analysis_numconv(w_ref, 7)

resol_list=["Eh192Ez30P3", "Eh384Ez60P3", "Eh1536Ez240P3"]
var_list={"Eh192Ez30P3": w, "Eh384Ez60P3": w_hi2, "Eh1536Ez240P3": w_hi4}
numerror_w = numconv.eval_numerror_conv(resol_list, var_list, w_ref, "l2")

numerror1 = numerror_w["Eh192Ez30P3"]
numerror2 = numerror_w["Eh1536Ez240P3"]
resol_fac = 8
conv_rate = np.log2(numerror1/numerror2)/np.log2(resol_fac)

os.makedirs(out_dir, exist_ok=True)
#---
fig, ax = plt.subplots(1,1,figsize=(9,5)) 
conv_rate_=conv_rate.rolling(center=True,x=7,z=1).mean()
w_ref_ = w_ref.copy()
conv_rate_.sel(z=slice(0e3,12e3),x=slice(95e3,145e3)).plot(vmin=1.5,vmax=4.5,ax=ax,cmap="jet")
cntr = w_ref_.sel(z=slice(0e3,12e3),x=slice(95e3,145e3)).plot.contour(ax=ax,levels=[-0.1,-0.075,-0.05,-0.025,0,0.025,0.05,0.075,0.1],colors='white')
cntr.clabel(fmt='%3i',fontsize=12,colors='white')
set_fig_XZ_axis(ax, True)
plt.savefig(f"{out_dir}/dist_conv_rate.png")

#---
fig, ax = plt.subplots(1,1,figsize=(9,5)) 
numerror_w_ = numerror_w["Eh384Ez60P3"].copy()
numerror_w_.rolling(center=True,x=1,z=1).mean().sel(z=slice(0e3,12e3),x=slice(95e3,145e3)).plot(ax=ax,cmap="jet",vmin=0,vmax=1e-4)
cntr = w_ref_.sel(z=slice(0e3,12e3),x=slice(95e3,145e3)).plot.contour(ax=ax,levels=[-0.1,-0.075,-0.05,-0.025,0,0.025,0.05,0.075,0.1],colors='white')
set_fig_XZ_axis(ax, True)
plt.savefig(f"{out_dir}/dist_L2error.png")
