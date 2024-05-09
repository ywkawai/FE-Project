import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import xarray as xr
import numpy as np
import os
import joblib

sys.path.append(os.path.join(os.path.dirname(__file__), '.'))
import analysis_common as ana_comm

class analysis_numconv:
  def __init__(self, dir, varname, postfix, pord=7) -> None:
    var_ref_tmp = ana_comm.open_merge_nc(dir, varname, postfix)
    self.node_num = pord+1
    if pord==7:
        w1d = [0.0357,0.210,0.341,0.412,0.412,0.341,0.21,0.0357]
    elif pord==3:
        w1d = [0.16666,0.83333,0.83333,0.16666]
        
    w_h1d = np.array(w1d*int(var_ref_tmp.shape[1]/self.node_num))
    w_v1d = np.array(w1d*int(var_ref_tmp.shape[0]/self.node_num))
    print(var_ref_tmp.shape)
    w = np.tile(w_h1d, (var_ref_tmp.shape[0],1))
    self.w_3D = w * np.transpose( np.tile(w_v1d, (var_ref_tmp.shape[1],1)) )
    self.fac=0.5*self.node_num

  def get_var_list( self, dir, resol_list, varname, time, isely=0):
      var_list = {}
      for resol in resol_list:
          var_list[resol] = ana_comm.open_merge_nc(f"{dir}/{resol}", varname, f"_y{isely}_t{time}") 
      return var_list

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

  def mkgraph_numconv( self, numerror1, numerror2, resol_fac, mask_fac=0.01, zlim=[0,30e3], xlim=[0e3,240e3], ax=None, vmax=4.5 ):
    if ax == None:
        f, ax = plt.subplots(1, 1,figsize=(6, 4))
    
    tmp = np.log2(numerror1/numerror2)/np.log2(resol_fac)
    maxval = np.abs(numerror2.sel(x=slice(xlim[0],xlim[1]),z=slice(zlim[0],zlim[1]))).max().values*mask_fac
    tmp = tmp.where(np.abs(numerror2) > maxval)
    tmp.rolling(center=True,x=16,z=1).mean().plot(ax=ax,vmin=1,vmax=vmax,cmap="jet",ylim=zlim,xlim=xlim)

  def print_numconv_info( self, numerror_list, xlim=[0,240e3], zlim=[0,15e3]):
    numerror_keys = list(numerror_list.keys())
    numerror_list = list(numerror_list.values())
    list_len = len(numerror_list)
    
    numconv_list = {}
    numerror_lc_list = {}    
    for i in range(0,list_len):
        numerror_lc_list[numerror_keys[i]] = numerror_list[i].sel(x=slice(xlim[0],xlim[1]),z=slice(zlim[0],zlim[1])).mean().values
        if i > 0:
            conv_rate = np.log2(numerror_list[i-1].sel(x=slice(xlim[0],xlim[1]),z=slice(zlim[0],zlim[1])).mean().values / numerror_list[i].sel(x=slice(xlim[0],xlim[1]),z=slice(zlim[0],zlim[1])).mean().values)
#            print(f"{numerror_keys[i]}: {l1error:.4E}   {conv_rate:.3f}")
        else:
            conv_rate = np.NaN
#            print(f"{numerror_keys[i]}: {l1error:.4E}   {conv_rate}")
        numconv_list[numerror_keys[i]] = conv_rate
        
    if list_len > 1:
        conv_rate = np.log2(numerror_list[0].sel(x=slice(xlim[0],xlim[1]),z=slice(zlim[0],zlim[1])).mean().values / numerror_list[list_len-1].sel(x=slice(xlim[0],xlim[1]),z=slice(zlim[0],zlim[1])).mean().values) / np.log2(2**(list_len-1))
        numconv_list["--"] = conv_rate

    for key, conv_rate in numconv_list.items():
        if key=="--":
            print(f"{key}: --  {numconv_list[key]:.3f}")                    
        else:
            print(f"{key}: {numerror_lc_list[key]:.4E}  {numconv_list[key]:.3f}")
    return numerror_lc_list, numconv_list

  def investigate_numconv( self, resol_list, var_listlist, var_ref_list, error_type, 
                          numconv_xlim, numconv_zlim, 
                          mkgraph_resol1, mkgraph_resol2, resol_fac, mask_fac, xlim, zlim, vmax ):

    var_num = len(var_listlist.values()) 
    print(var_num)   
    f, axes = plt.subplots(var_num,1,figsize=(6, 4*var_num))
    for i, (varname, var_list) in enumerate(var_listlist.items()):
        numerror_list = self.eval_numerror_conv( resol_list, var_list, var_ref_list[varname], error_type )
        self.mkgraph_numconv( numerror_list[mkgraph_resol1], numerror_list[mkgraph_resol2], resol_fac, mask_fac, xlim=xlim, zlim=zlim,
                        ax=axes[i], vmax=vmax )    
        print(f"* {varname} =============")
        numerror_lc, numconv = self.print_numconv_info( numerror_list, xlim=numconv_xlim, zlim=numconv_zlim)
        
  def output_numconv_data( self, resol_list, varname_list, error_type,
                          numconv_xlim, numconv_zlim, time_list, tmp_dir, tmp_dir_ref ):
    resol_list2 = [*resol_list, "--"]

    ds = xr.Dataset()
    for varname in varname_list:        
        np_numerror_lc = np.zeros((len(time_list),len(resol_list)))
        np_numconv_lc = np.zeros((len(time_list),len(resol_list2)))
        
        for t, time in enumerate(time_list): 
            var_list = self.get_var_list( tmp_dir, resol_list, varname, time )
            var_ref = ana_comm.open_merge_nc( tmp_dir_ref, varname, f"_y0_t{time}")
            numerror_list = self.eval_numerror_conv( resol_list, var_list, var_ref, error_type )
            numerror_lc, numconv = self.print_numconv_info( numerror_list, xlim=numconv_xlim, zlim=numconv_zlim)
            np_numerror_lc[t,:] = list(numerror_lc.values())
            np_numconv_lc[t,:] = list(numconv.values())
            
        xr_numerror_lc = xr.DataArray( np_numerror_lc, coords={'time': time_list, 'resol': resol_list}, dims=('time', 'resol') ).rename(f"{varname}_{error_type}")
        xr_numconv_lc = xr.DataArray( np_numconv_lc, coords={'time': time_list, 'resol2': resol_list2}, dims=('time', 'resol2') ).rename(f"{varname}_{error_type}_convrate")        
        ds[f"{varname}_{error_type}"] = xr_numerror_lc
        ds[f"{varname}_{error_type}_convrate"] = xr_numconv_lc

    ds.to_netcdf(f"{tmp_dir}/numerror_{error_type}.nc") 