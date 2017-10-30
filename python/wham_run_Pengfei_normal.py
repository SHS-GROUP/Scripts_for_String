#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
from __future__ import print_function
from wham_v1_beta5 import *
from matplotlib import rc
import matplotlib.pyplot as plt
from optparse import OptionParser

parser = OptionParser("wham_run_Pengfei_normal.py -i image_num -d dim_num -c cycle_num -b bin_num --bs bin_size --bt btstrap --fx fx_file --dir dirpath")
parser.set_defaults(binsize=0.0, binnum=0, btstrap=0, dirpath='.', fxfile='')
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
parser.add_option("--bt", dest="btstrap", type='int',
                  help="bootstrap is used (1) or not (0) [default: 0]")
parser.add_option("--fx", dest="fxfile", type='string',
                  help="The WHAM weights file")
parser.add_option("--dir", dest="dirpath", type='string',
                  help="The directory path of the data file")

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

if options.btstrap not in [0, 1]:
    raise ValueError('btstrap need to be 0 or 1!')
btstrap = options.btstrap #Whether to use fake boots trapping

react_paths = [r'CO ($\AA$)', r'OH ($\AA$)', r'CH ($\AA$)'] #Reaction path names
wham_conv = 0.001 #WHAM converge creteria

if options.dirpath[-1] != '/':
    dirpath = options.dirpath + '/'

#
# Do the WHAM cycle
#
if options.fxfile == '':
    data_dict, num_sims, string_seq, expFx, Ubiasl = wham(dirpath, dim, react_paths, first_num_imgs, num_cycles, wham_conv, btstrap)
else:
    data_dict, num_sims, string_seq, expFx, Ubiasl = wham(dirpath, dim, react_paths, first_num_imgs, num_cycles, wham_conv, btstrap, options.fxfile)

#
# For each dimension
#
for i in xrange(0, dim):
    plot_dim = i + 1
    figname = 'R' + str(i+1) + '.pdf'
    gene_free_ene_1D(data_dict, num_sims, dim, expFx, Ubiasl, [plot_dim], [1.0], num_bins, figname, react_paths[i])

#
# For linear combination of 1D system: e.g. R3-R2
#
gene_free_ene_1D(data_dict, num_sims, dim, expFx, Ubiasl, [3, 2], [1.0, -1.0], num_bins, 'R3-R2.pdf', r'CH-OH ($\AA$)')

#
# For normal 2D system: e.g. R2 vs R1
#
gene_free_ene_2D(data_dict, num_sims, dim, expFx, Ubiasl, [2], [1.0], [1], [1.0], num_bins, num_bins, 'R2_R1.pdf', r'OH ($\AA$)', r'CO ($\AA$)', string_seq, 'rainbow')

#
# For normal 2D system: e.g. R3 vs R1
#
gene_free_ene_2D(data_dict, num_sims, dim, expFx, Ubiasl, [3], [1.0], [1], [1.0], num_bins, num_bins, 'R3_R1.pdf', r'CH ($\AA$)', r'CO ($\AA$)', string_seq, 'rainbow')

#
# For 2D system with linear combination: e.g. R3-R2 vs R1
#
gene_free_ene_2D(data_dict, num_sims, dim, expFx, Ubiasl, [3, 2], [1.0, -1.0], [1], [1.0], num_bins, num_bins, 'R3-R2_R1.pdf', r'CH-OH ($\AA$)', r'CO ($\AA$)', string_seq, 'rainbow')

#
# For 2D system with linear combination: e.g. R1-R2-R3 vs R1
#
#gene_free_ene_2D(data_dict, num_sims, dim, expFx, Ubiasl, [1, 2, 3], [1.0, -1.0, -1.0], [1], [1.0], num_binsX, num_binsY, 'R1-R2-R3_R3.pdf', r'CO-OH-CH ($\AA$)', r'CO ($\AA$)', string_seq, 'rainbow')

cost_time = time.time() - start_time0
print("It costs %f seconds to finish the job!" %cost_time)
quit()

