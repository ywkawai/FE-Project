import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

def get_moni_data(dir, dtsec, t_offset=0.0):
    skipfooter=3
    #names=['dummy1', 'step', 'dummy2', 'DDENS', 'QV', 'QC', 'QR', 'QTOT', 'EVAP', 'PREC', 'ENGT', 'ENGP', 'ENGK', 'ENGI' ]
    names=['dummy1', 'step', 'dummy2', 'DDENS', 'ENGT', 'ENGP', 'ENGK', 'ENGI' ]

    df = pd.read_csv(f'{dir}/monitor.peall', names=names, skiprows=range(0, 1), delim_whitespace=True)
    df = df.drop(columns=['dummy1', 'dummy2'])
    df.loc[:,"sec"] = t_offset + (df.loc[:,"step"]-1)*dtsec
    df = df.set_index('sec')
    ds = xr.Dataset.from_dataframe(df)
    return ds

def get_moni_data_list(exp_dir, run_num, dtsec, run_num0=1, t_offset0=0.0, val_offset0=0.0):
    ds_list = []
    dir_list = []
    for runno in range(0,run_num):
        dir_list.append(f"{exp_dir}/run{run_num0+runno}")
    
    i = 1
    t_offset = t_offset0
    val_offset = val_offset0
    
    for dir in dir_list:
        ds_ori = get_moni_data(dir, dtsec, t_offset)
        time = ds_ori["sec"].values
        tlen = len(time)    
        if i > 1:
            ds = ds_ori.isel(sec=slice(1,tlen)) + val_offset
        else:
            ds = ds_ori + val_offset0
        
        ds_list.append(ds)
        i = i+1
        t_offset = time[tlen-1]
        val_offset = val_offset + ds_ori.isel(sec=tlen-1)

    return xr.concat(ds_list, "sec")
    
def mkgraph(ds_list, varname, ylim, color_list, linestyle_list):
    fig, ax = plt.subplots(figsize=(8,3))
    ax.set_ylim(ylim)
    ax.set_ylabel(varname)
    ax.set_xlabel("time [sec]")
    #ax.set_xlim([0,3e4])
    for key, ds in ds_list.items():
        time = ds.sec
        ax.plot(time, ds[varname], label=key, color=color_list[key], linestyle=linestyle_list[key])
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
