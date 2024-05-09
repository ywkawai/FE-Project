import xarray as xr
import numpy as np
import mkgraph_common as mkgraph
import matplotlib.ticker as tick
import matplotlib.pyplot as plt

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
  panelID,PRC_NUM_X,PRC_NUM_Y, 
  PRC_NUM_X_s,PRC_NUM_X_e,PRC_NUM_Y_s,PRC_NUM_Y_e):
  prc_offset = (panelID-1)*PRC_NUM_X*PRC_NUM_Y  
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

def set_fig_LonZ_axis(ax, xlabel):
  ax.tick_params(labelsize=18, length=10)  
#  print("ylabel:", xlabel)
  ax.xaxis.set_major_locator(tick.MultipleLocator(30.0))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(10.0))    
  ax.xaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis0_fmt))
  if (xlabel):  
    ax.set_xlabel('longitude [deg]', fontsize=18)  
  ax.yaxis.set_major_locator(tick.MultipleLocator(2000.0))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(500.0))  
  ax.yaxis.set_major_formatter(tick.FuncFormatter(vert_1Daxis_fmt))
  ax.set_ylabel('height [km]', fontsize=18)

#------------------

fpath = get_fpathlist(target_dir,"Eh24Ez12P3/outdata_analysis/history","","",3,16,16, 0,16,7,8)
df = merge_hist(fpath)
fpath_hi = get_fpathlist(target_dir,"Eh96Ez36P3/outdata_analysis/history","","",3,16,16, 0,16,7,8)
df_hi = merge_hist(fpath_hi)

fpath_ref = get_fpathlist(target_dir,"Eh96Ez36P7/history","","",3,16,16, 0,16,7,8)
df_ref = merge_hist(fpath_ref)

w = df["W"].isel(y=47,time=12)
w_hi = df_hi["W"].isel(y=47,time=12)
w_ref = df_ref["W"].isel(y=47,time=12)
numconv = analysis_numconv(w_ref, 7)

resol_list=["Eh24Ez12P3", "Eh96Ez36P3"]
var_list={"Eh24Ez12P3": w, "Eh96Ez36P3": w_hi}
numerror_w = numconv.eval_numerror_conv(resol_list, var_list, w_ref, "l2")

numerror1 = numerror_w["Eh24Ez12P3"]
numerror2 = numerror_w["Eh96Ez36P3"]
resol_fac = 4
conv_rate = np.log2(numerror1/numerror2)/np.log2(resol_fac)
Radius=38219.6760

#---
fig, ax = plt.subplots(1,1,figsize=(9,5)) 
conv_rate_=conv_rate.rolling(center=True,x=7,z=1).mean()
conv_rate_.coords["x"] = conv_rate_.coords["x"]*180/np.pi + 180
w_ref_ = w_ref.copy()
w_ref_.coords["x"] = w_ref_.coords["x"]*180/np.pi + 180
conv_rate_.sel(z=slice(60e3,72e3),x=slice(140,220)).plot(vmin=1.5,vmax=4.5,ax=ax,cmap="jet")
cntr = w_ref_.sel(z=slice(60e3,72e3),x=slice(140,220)).plot.contour(ax=ax,levels=[-0.1,-0.075,-0.05,-0.025,0,0.025,0.05,0.075,0.1],colors='white')
set_fig_LonZ_axis(ax, True)
plt.savefig(f"{out_dir}/dist_conv_rate.png")

#---
fig, ax = plt.subplots(1,1,figsize=(9,5)) 
numerror_w_ = numerror_w["Eh96Ez36P3"].copy()
numerror_w_.coords["x"] = numerror_w_.coords["x"]*180/np.pi + 180
numerror_w_.rolling(center=True,x=1,z=1).mean().sel(z=slice(60e3,72e3),x=slice(140,220)).plot(ax=ax,cmap="jet",vmin=0,vmax=1e-4)
cntr = w_ref_.sel(z=slice(60e3,72e3),x=slice(140,220)).plot.contour(ax=ax,levels=[-0.1,-0.075,-0.05,-0.025,0,0.025,0.05,0.075,0.1],colors='white')
set_fig_LonZ_axis(ax, True)
plt.savefig(f"{out_dir}/dist_L2error.png")
