#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
from __future__ import print_function
from numpy import array, mean, std
from optparse import OptionParser
from matplotlib import pyplot as plt
from matplotlib import rc
from string_functions import read_list, read_2d_free_energy, write_xy_lists, write_2d_free_energy

rc("font", **{"family": 'sans-serif', "sans-serif": ['Arial'], "size": 20, "weight": 'normal'})

parser = OptionParser("Usage: cal_string_avg.py -i input_name_lists -o output_file -d dimension")
parser.add_option("-i", dest="inputf", type='string',
                  help="Input file containing the name lists")
parser.add_option("-o", dest="outputf", type='string',
                  help="Output file name")
parser.add_option("-d", dest="dim", type='int',
                  help="Topology file name")
(options, args) = parser.parse_args()

name_list = []
r_inputf = open(options.inputf, 'r')
for rline in r_inputf:
    line = rline.strip('\n')
    line = line.strip(' ')
    name_list.append(line)
r_inputf.close()

if options.dim == 1:
    data_super_list = []
    for name in name_list:
        data_list = read_list(name, 1)
        data_super_list.append(data_list)

    data_super_list = array(data_super_list)
    avg_data_list = mean(data_super_list, axis=0)
    std_data_list = std(data_super_list, axis=0)

    # Print out the free energy profile and errorbars
    write_xy_lists(options.outputf, avg_data_list, std_data_list)

    # Plot the free energy profile with error bar
    x_list = xrange(1, len(data_list)+1)
    plt.errorbar(x_list, avg_data_list, yerr=std_data_list, linewidth=2.0)
    plt.axes().xaxis.get_major_locator().set_params(nbins=5)
    plt.axes().yaxis.get_major_locator().set_params(nbins=5)
    plt.axes().tick_params(axis='x', which='minor', length=2)
    plt.ylabel('Free Energy (kcal/mol)', weight='bold', fontsize=15)
    plt.xlabel('Image Number', weight='bold', fontsize=15)
    plt.tick_params(which='both', width=2)
    plt.tight_layout()
    plt.savefig(options.outputf+'.pdf')
    plt.close()

elif options.dim == 2:
    data_super_list = []
    for name in name_list:
        data_list = read_2d_free_energy(name)
        data_super_list.append(data_list)

    data_super_list = array(data_super_list)
    avg_data_list = mean(data_super_list, axis=0)
    std_data_list = std(data_super_list, axis=0)

    # Print out the 2d free energy with errorbars
    write_2d_free_energy(options.outputf, avg_data_list)
    write_2d_free_energy(options.outputf+'.errbar', std_data_list)

    # There is no error bar plotting for 2D free energy surface

quit()

