import os
import sys

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

#------------------------------
dir = "../Dx3.1m_P7/run51/outdata/"
PRC_NUM_X = 16; PRC_NUM_Y = 8;

out_pngdir = "./fig_snapshot/"

#------------------------------
def get_fpathlist(dir1,ftype,domainlabel,timelabel,
                  PRC_NUM_X,PRC_NUM_Y, PRC_NUM_XS,PRC_NUM_XE,PRC_NUM_YS,PRC_NUM_YE):
    return [ [ dir1 + "{0}{1}{2}.pe{3:06d}.nc".format( ftype, domainlabel, timelabel, i + j*PRC_NUM_X ) for i in range(PRC_NUM_XS-1,PRC_NUM_XE) ] for j in range(PRC_NUM_YS-1,PRC_NUM_YE) ]

def merge_xy(fpath, dim=["y","x"]):
    return xr.open_mfdataset(fpath, decode_times=False,combine="nested", concat_dim=dim)

def open_tmp_nc(dir, varname):
  print(f'{dir}/history.pe*.nc')
  return xr.open_mfdataset(f'{dir}/history.pe*.nc', decode_times=False, combine='by_coords')[varname]
#------------------------------

fpath = get_fpathlist(dir, "history", "","", PRC_NUM_X,PRC_NUM_Y, 1,PRC_NUM_X,1,PRC_NUM_Y)
ds = merge_xy(fpath)

# --- plot ---
if not os.path.exists(out_pngdir):
    os.makedirs(out_pngdir)


fig, ax = plt.subplots(1,1,figsize=(7.5,6))
ds.isel(time=0).interp(z=800)["W"].plot(cmap="RdBu_r", vmin=-6,vmax=6, ax=ax)
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
ax.set_title("")
plt.savefig(f"{out_pngdir}/W_snapshot.png", bbox_inches='tight')

fig, ax = plt.subplots(1,1,figsize=(7,3))
(ds.sel(time=0).interp(y=800)["U"]-5.0).plot(cmap="RdBu_r", vmin=-6,vmax=6, ax=ax)
ax.set_xlabel("y [m]")
ax.set_ylabel("z [m]")
ax.set_title("")
plt.savefig(f"{out_pngdir}/U_yz_snapshot.png", bbox_inches='tight')

fig, ax = plt.subplots(1,1,figsize=(7,3))
ds.sel(time=0).interp(y=800)["W"].plot(cmap="RdBu_r", vmin=-6,vmax=6, ax=ax)
ax.set_xlabel("y [m]")
ax.set_ylabel("z [m]")
ax.set_title("")
plt.savefig(f"{out_pngdir}/W_yz_snapshot.png", bbox_inches='tight')

fig, ax = plt.subplots(1,1,figsize=(7,3))
ds.sel(time=0).interp(y=800)["V"].plot(cmap="RdBu_r", vmin=-6,vmax=6, ax=ax)
ax.set_xlabel("y [m]")
ax.set_ylabel("z [m]")
ax.set_title("")
plt.savefig(f"{out_pngdir}/V_yz_snapshot.png", bbox_inches='tight')

fig, ax = plt.subplots(1,1,figsize=(7,3))
ds.isel(time=0).interp(y=800)["PT"].plot(cmap="jet", vmin=299.8,vmax=300.9, ax=ax)
ax.set_xlabel("y [m]")
ax.set_ylabel("z [m]")
ax.set_title("")
plt.savefig(f"{out_pngdir}/PT_yz_snapshot.png", bbox_inches='tight')