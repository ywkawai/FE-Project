import sys
import analysis_monitor_file as analysis_moni
import argparse

parser = argparse.ArgumentParser(description='analysis_monitor.')
parser.add_argument('inputFilePath', help='Input file path')
parser.add_argument('varName', help="Variable name")
parser.add_argument('DelTimeSec', help="Timestep [sec] of model", type=float)
parser.add_argument('outputFilePath', help='Output file path (figure)')
parser.add_argument('--ylim', help='Set ylim in figure', nargs=2, type=float, default=[0.0,0.0] )

args = parser.parse_args()
moni_file = args.inputFilePath
dtsec     = args.DelTimeSec
vname     = args.varName
fig_file  = args.outputFilePath
ylim      = args.ylim

print(f"Target file: {moni_file}")
print(f"DelTimeSec: {dtsec}")
print(f"ylim: {ylim}")

ds = analysis_moni.get_moni_data(moni_file, dtsec)
da_list = {}; color_list = {}; linestyle_list = {}
da_list[vname] = ds[vname]
color_list[vname] = "black"
linestyle_list[vname] = "-"

if (ylim[0]==ylim[1]):
    minv = ds[vname].values.min()
    maxv = ds[vname].values.max()
    ylim[0] = minv - (maxv-minv)*0.02
    ylim[1] = maxv + (maxv-minv)*0.02

fig = analysis_moni.mkgraph( da_list, vname, ylim, color_list, linestyle_list)
print(f"Output figure: {fig_file}")
fig.savefig(fig_file, bbox_inches='tight')