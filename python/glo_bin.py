#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# File name glo_bin.py
from string_functions import get_color_dict, read_dat_file
from numpy import linspace, histogram
import matplotlib.pyplot as plt
from optparse import OptionParser

parser = OptionParser("Usage: glo_bin.py -c cycle -w num_windows -b num_bins -d dim")
parser.add_option("-c", dest="cycle", type='int', help="Cycle number")
parser.add_option("-w", dest="windows", type='int', help="Number of windows in the cycle")
parser.add_option("-b", dest="bins", type='int', help="Number of bins used")
parser.add_option("-d", dest='dim', type='int', help="The dimension to check")
(options, args) = parser.parse_args()

datal = []
for i in range(1, options.windows+1):
    datf = './optstring_calcs' + str(options.cycle) + '/s' + str(options.cycle) + '_i' + str(i) + '_collect.dat'
    data = read_dat_file(datf)
    datal.append(data)

# Get the sampling data for a certain dimension
data_RC = []
data_RCl = []
for i in xrange(0, len(datal)):
    data_RC_per_sim = []
    for j in xrange(0, len(datal[i])):
        data_val = datal[i][j][options.dim-1]
        data_RC.append(data_val)
        data_RC_per_sim.append(data_val)
    data_RCl.append(data_RC_per_sim)

# Determine the plotting bin size
val_min = min(data_RC) - 0.00001
val_max = max(data_RC) + 0.00001

bins_RC = linspace(val_min, val_max, options.bins+1)
binsize_RC = bins_RC[1] - bins_RC[0]
plot_bins_RC = bins_RC[0:-1]
xaxis_RC = [i + binsize_RC/2.0 for i in plot_bins_RC]

# Global binning
color_dict = get_color_dict(len(datal))

for i in xrange(0, len(datal)):
    wlabel = 'w' + str(i+1)
    data_RC_per_sim = data_RCl[i]
    count_RC, bins_RCp = histogram(data_RC_per_sim, bins=bins_RC)
    clr = color_dict[i+1]
    plt.plot(xaxis_RC, count_RC, color=clr, linestyle='-', linewidth=2.0, label=wlabel)
    plt.legend(loc=1)

figname = 'cycle' + str(options.cycle) + '_' + 'R' + str(options.dim) + '_globin.pdf'
plt.savefig(figname)
plt.close()

quit()

