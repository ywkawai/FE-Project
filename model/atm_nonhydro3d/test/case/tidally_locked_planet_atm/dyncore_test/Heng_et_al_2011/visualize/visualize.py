import numpy as np
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
from joblib import Parallel, delayed

TMP_DIR_SUFFIX="./tmp_outdata"
PRC_NUM_X=4; PRC_NUM_Y=4
OUT_FIG_dir=f"./"
Nproc = 4

def get_fpathlist(dir1,ftype,domainlabel,timelabel,PRC_NUM_X,PRC_NUM_Y):
  return [ [ dir1 + "{0}{1}{2}.pe{3:06d}.nc".format( ftype, domainlabel, timelabel, i + j*PRC_NUM_X ) for i in range(PRC_NUM_X) ] for j in range(PRC_NUM_Y) ]

def merge_xy(fpath, dim=["y","x"]):
  return xr.open_mfdataset(fpath, decode_times=False,combine="nested", concat_dim=dim)

def output_tmean(ds, plev_hPa_list, varname_list, outdir):
    r = Parallel(n_jobs=Nproc)( [delayed(output_tmean_sub)(ds, plev, varname, outdir) for plev in plev_hPa_list for varname in varname_list] )

def output_tmean_sub(ds, pres_hPa, varname, outdir):
  print(f"output_tmean: plev={pres_hPa}hPa varname={varname} ..")
  ds_tmp = ds[varname].sel(p=pres_hPa*1e2)
  if ds_tmp.p.size > 1:
    ds_tmp = ds_tmp.isel(p=0)
  ds_tmp.mean(["time"]).to_netcdf(f"{outdir}/{varname}_{pres_hPa}hPa.nc")
    
def get_var_tmean(varname, runno_list, plev_hPa):
  var = xr.open_mfdataset(f"{TMP_DIR_SUFFIX}1/{varname}_{plev_hPa}hPa.nc", decode_times=False)[varname]
  for i in range(1,len(runno_list)):
    runno = runno_list[i]
    var = var + xr.open_mfdataset(f"{TMP_DIR_SUFFIX}{runno}/{varname}_{plev_hPa}hPa.nc", decode_times=False)[varname]    
  return var / float(len(runno_list))

def mkgraph_uvt_horidist(plev, runno_list, pngname, levels):
  umet = get_var_tmean("Umet", runno_list, plev)
  vmet = get_var_tmean("Vmet", runno_list, plev)
  temp = get_var_tmean("T", runno_list, plev)
  
  skip=16
  fig, ax = plt.subplots(1,1,figsize=(10,5))
  ax.set_ylim(-87.5,87.5)
  temp.plot.contourf(ax=ax,cmap="jet", levels=levels)
  ax.quiver(temp.lon[::skip],temp.lat[::skip],umet[::skip,::skip],vmet[::skip,::skip])  
  plt.savefig(pngname)

#-----
plev_hPa_list = [250, 680, 950, 1000]
varname_list = ["Umet", "Vmet", "T"]
temp_level_list = {
  250: np.arange(225,230,0.5),
  680: np.arange(290,304,1), 
  950: np.arange(265,335,5), 
  1000: np.arange(255,335,5), 
}
runno_list = [1 ,2]

for runno in runno_list:  
  ds = merge_xy( get_fpathlist(f"../run{runno}/outdata_p/", "history","","", 4,4), dim=["lat", "lon"])

  out_dir = f"{TMP_DIR_SUFFIX}{runno}"
  os.makedirs(out_dir, exist_ok=True)
  print(f"output_tmean: runno={runno}")  
  output_tmean(ds, plev_hPa_list, varname_list, out_dir)
  
for plev_hPa in plev_hPa_list:
  print(f"mkgraph: plev={plev_hPa}hPa")
  mkgraph_uvt_horidist(plev_hPa, runno_list, 
                       f"{OUT_FIG_dir}/uvt_{plev_hPa}hPa.png", temp_level_list[plev_hPa])

  
