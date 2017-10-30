#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
from __future__ import print_function
from wham_v1_beta5 import *
from matplotlib import rc
import matplotlib.pyplot as plt
from optparse import OptionParser

parser = OptionParser("wham_run_Pengfei_avg.py -i image_num -d dim_num -c cycle_num -b bin_num --bs bin_size --fx fx_file --f2avg file2avg --fclm file_column -o outputf --dir dirpath")
parser.set_defaults(binsize=0.0, binnum=0, dirpath='.', fxfile='', fclm=1)
parser.add_option("-i", dest="image", type='int',
                  help="Image number of the initial string")
parser.add_option("-d", dest="dim", type='int',
                  help="Dimension number")
parser.add_option("-c", dest="cycle", type='int',
                  help="Cycle number")
parser.add_option("-b", dest="binnum", type='int',
                  help="bin number (need to specify either binnum or binsize)")
parser.add_option("--bs", dest="binsize", type='float',
                  help="bin size (need to specify either binnum or binsize)")
parser.add_option("--fx", dest="fxfile", type='string',
                  help="The WHAM weights file [default: none, will calculate automatically, but it will cost a lot of time!]")
parser.add_option("--f2avg", dest="f2avg", type='string',
                  help="Data file for averaging")
parser.add_option("--fclm", dest="fclm", type='int',
                  help="Column file of the data file to analysis [default: 1]")
parser.add_option("-o", dest="output", type='string',
                  help="Output file name")
parser.add_option("--dir", dest="dirpath", type='string',
                  help="The directory path of the data file [default: current directory]")
(options, args) = parser.parse_args()

###############################################################################
                                   #MAIN PROGRAM
###############################################################################
start_time0 = time.time() #Count the time

#
# Setting variables
#
dim=options.dim #Dimension of the reaction coordinates
first_num_imgs=options.binnum #Number of images for the string
num_cycles=options.cycle #Number of iteration for the string calculations

if options.binnum == 0 and options.binsize == 0.0:
    raise ValueError('Need to specify either binsize or binnum!')
elif options.binnum != 0 and options.binsize != 0.0:
    raise ValueError('Need to specify either binsize or binnum, not both of them!')
elif options.binsize == 0.0:
    num_bins=[options.binnum, False] #Number of bins, specify the bin number
else:
    num_bins = [False, options.binsize] # Specify the bin size

react_paths = [r'CO ($\AA$)', r'OH ($\AA$)', r'CH ($\AA$)'] #Reaction path names
wham_conv = 0.001 #WHAM converge creteria

if options.dirpath[-1] != '/':
    dirpath = options.dirpath + '/'

#
# Do the WHAM cycle
#
btstrap = 0
if options.fxfile == '':
    data_dict, num_sims, string_seq, expFx, Ubiasl = wham(dirpath, dim, react_paths, first_num_imgs, num_cycles, wham_conv, btstrap)
else:
    data_dict, num_sims, string_seq, expFx, Ubiasl = wham(dirpath, dim, react_paths, first_num_imgs, num_cycles, wham_conv, btstrap, options.fxfile)

data2avg = read_list(options.f2avg, options.fclm)
gene_free_ene_2D(data_dict, num_sims, dim, expFx, Ubiasl, [3], [1.0], [1], [1.0], [False, 0.1], [False, 0.1], options.output, react_paths[2], react_paths[0], string_seq, 'rainbow', 1, data2avg)

cost_time = time.time() - start_time0
print("It costs %f seconds to finish the job!" %cost_time)
quit()

