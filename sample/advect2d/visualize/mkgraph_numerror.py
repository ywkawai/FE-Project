import sys
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr

args = sys.argv
in_numerror_logfile = args[1]
outpng_basename     = args[2]

def mkgraph_numerror( time, numerror, fig_title, ylbl, ylim_min, ylim_max, out_fname ):
  fig, ax = plt.subplots(figsize=(6,4))

#  ax.set_xlim(0.0, 1.0)
  ax.plot( time, numerror )
  ax.set_ylim(ylim_min , ylim_max)
  ax.set_yscale('log')
  ax.set_xlabel("time [s]")
  ax.set_ylabel(ylbl)
  ax.grid()
  ax.set_title(fig_title) 
  plt.savefig( out_fname )

names=['dummy1', 'tsec', 'L1_error', 'L2_error', 'Linf_error', 'Ediss', 'Edisp']

df = pd.read_csv(in_numerror_logfile, names=names, skiprows=range(0, 1), sep='\s+')
df = df.drop(columns=['dummy1']).set_index('tsec')
ds = xr.Dataset.from_dataframe(df)

time = ds["tsec"]
mkgraph_numerror( time, ds["L1_error"], "L1 error", "L1 error", 1e-9, 1e-2, f"{outpng_basename}_l1.png")
mkgraph_numerror( time, ds["L2_error"], 'L2 error', "L2 error", 1e-9, 1e-2, f"{outpng_basename}_l2.png")
mkgraph_numerror( time, ds["Linf_error"], 'Linf error', "Linf error", 1e-8, 1e-1, f"{outpng_basename}_linf.png")
