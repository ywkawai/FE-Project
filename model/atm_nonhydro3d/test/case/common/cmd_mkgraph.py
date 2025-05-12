import sys
import mkgraph_plot as mkgraph
import argparse
import re
import xarray as xr
import numpy as np

def parse_gturl(gturl):
    pattern = r'^(.*)\.(pe[0-9\*]+)\.nc@([^,]+)(?:,(.*))?$'
    match = re.match(pattern, gturl)

    if match:
        filename = match.group(1)       # "history.pe0.nc"
        peinfo   = match.group(2)
        variable = match.group(3)       # "U"
        params_str = match.group(4)     # "x=0,t=0" or None

        params = {}
        if params_str:
            param_pairs = [p.split('=') for p in params_str.split(',') if '=' in p]
            for k, v in param_pairs:
                try:
                    params[k] = float(v)
                except ValueError:
                    params[k] = v
        return filename, variable, params
    else:
        print("Error! Check arguments")
        exit(1)
    
def get_fpathlist(fname,
                  PRC_NUM_X,PRC_NUM_Y, PRC_NUM_XS=1,PRC_NUM_XE=-1,PRC_NUM_YS=1,PRC_NUM_YE=-1):
    if PRC_NUM_XE == -1:
        PRC_NUM_XE = PRC_NUM_X    
    if PRC_NUM_YE == -1:
        PRC_NUM_YE = PRC_NUM_Y
    return [ [ "{0}.pe{1:06d}.nc".format( fname, i + j*PRC_NUM_X ) for i in range(PRC_NUM_XS-1,PRC_NUM_XE) ] for j in range(PRC_NUM_YS-1,PRC_NUM_YE) ]

def merge_xy(fpath, dim=["y","x"]):
    return xr.open_mfdataset(fpath, decode_times=False,combine="nested", concat_dim=dim)

def get_index(var, axis_name, val):
    ax = var.coords[axis_name]
    return np.abs(ax-val).argmin()

#----------------------
parser = argparse.ArgumentParser(description='mkgraph')
parser.add_argument('gturl', help='gturl')
parser.add_argument('outputFilePath', help='Output file path (figure)')
parser.add_argument('--range', help='Set range', nargs=2, type=float, default=[0.0,0.0] )
parser.add_argument('--figsize', help='Figure size', nargs=2, type=int, default=None )
parser.add_argument('--cmap', help='Set colormap', default="jet" )
parser.add_argument('--prc_num_xy', help='Set PRC_NUM_X, PRC_NUM_Y', nargs=2, type=int, default=[1,1] )
parser.add_argument('--title', help='Figure title', default=None )
parser.add_argument('--exch', help='Exchange horizontal and vertical axes', action='store_true' )

args = parser.parse_args()
filename, variable, params = parse_gturl(args.gturl)
output_fpath = args.outputFilePath

prc_num_x = args.prc_num_xy[0]
prc_num_y = args.prc_num_xy[1]

vmin = args.range[0]
vmax = args.range[1]
if (vmin == vmax):
    vmin = None; vmax = None
    
cmap = args.cmap

print("File name:", filename)
print("Variable name:", variable)
print("Parameters:", params)
fpath = get_fpathlist( filename, prc_num_x, prc_num_y )
da = merge_xy(fpath)[variable]

for key, v in params.items():
#    print(f"Sel: {key}={v}")
    ind = get_index(da, key, v)
    da = da.isel({key: ind})

if args.title:
    title = args.title
else:
    title = f"{variable} [{da.units}]"

fig = mkgraph.plot(da, vmin=vmin, vmax=vmax, cmap=cmap, title=title, figsize=args.figsize, exch=args.exch)
fig.savefig(output_fpath, bbox_inches='tight')