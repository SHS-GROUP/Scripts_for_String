#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
from __future__ import print_function
from math import sqrt
from string_functions import read_list
from numpy import log10
from matplotlib import pyplot as plt
from optparse import OptionParser
from kie_functions import *

parser = OptionParser("Usage: cal_kie.py -f R1_free_ene_file -c code_mode -p para_file [-w Wavefunction]")

parser.set_defaults(cmode='python', wfn='ms')

parser.add_option("-f", dest="r1fef", type='string',
                  help="R1 free energy file")
parser.add_option("-c", dest="cmode", type='string',
                  help="Code mode (math or python)")
parser.add_option("-p", dest="pfile", type='string',
                  help="Parameter file")
parser.add_option("-w", dest="wfn", type='string',
                  help="Wavefunction: ms, ho, or bs")
parser.add_option("--t1", dest="tfname1", type='string',
                  help="File containig 1000/T and Hrate")
parser.add_option("--t2", dest="tfname2", type='string',
                  help="File containig 1000/T and KIE")
(options, args) = parser.parse_args()

###############################################################################
#                                    Constants
###############################################################################

T_list1 = read_list(options.tfname1, 1)
kh_list = read_list(options.tfname1, 2)

T_list2 = read_list(options.tfname2, 1)
kd_list = read_list(options.tfname2, 2)
KIE_list = read_list(options.tfname2, 3)

#WT_HT_list = [3.5971, 3.5336, 3.4722, 3.4130, 3.3557, 3.3003, 3.2468, 3.1949, 3.1447, 3.0960]
#WT_Hrate_list = [213.35, 228.94, 269.96, 275.10, 327.19, 297.11, 329.06, 301.77, 348.08, 327.51]

#WT_HT_list = [3.5971, 3.5336, 3.4722, 3.3557, 3.3003, 3.2468, 3.1949, 3.1447]
#WT_Hrate_list = [213.35, 228.94, 269.96, 327.19, 297.11, 329.06, 301.77, 348.08]

#WT_DT_list = [3.5971, 3.5336, 3.4722, 3.3557, 3.3003, 3.2468, 3.1949, 3.1447]
#WT_Drate_list = [2.161, 2.4767, 2.9864, 4.2919, 3.6822, 5.023, 4.2442, 4.0017]
#WT_KIE_list = [WT_Hrate_list[i]/WT_Drate_list[i] for i in xrange(len(WT_Hrate_list))]

# The point with temperature as 1000/3.3 = 303 K
abs_T_list1 = [abs(i-1000.0/303.0) for i in T_list1]
temp2fit = abs_T_list1.index(min(abs_T_list1))

###############################################################################
#                                  Main program
###############################################################################
#
# Read the W(R) list
#
R_list = read_list(options.r1fef, 1)
#dR = R_list[1] - R_list[0]
WR_list = read_list(options.r1fef, 2)
para_list = read_para_file(options.pfile)

# para_list for Morse: R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el, lamb
# para_list for HO: R, freq1, req1, freq2, req2, dG0, V_el, lamb
# para_list for Bspline: R, dG0, lamb

umax = 3
vmax = 3
Qm1 = 1.0

# Get the parition function parameter
T = 1000.0 / T_list1[temp2fit]

if options.wfn.lower() == "ms":
    k_h, k_d = get_ks_Morse(options.cmode, kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1, 1)
elif options.wfn.lower() == "ho":
    k_h, k_d = get_ks_HO(kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1, 1)
elif options.wfn.lower() == "bs":
    npots = 128
    k_h, k_d = get_ks_bspline(kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, npots, Qm1, 1)

Qm1 = kh_list[temp2fit] / k_h

"""
# Get the lowest temperature
T = 1000.0 / WT_HT_list[0]

if options.wfn.lower() == "ms":
    k_h, k_d = get_ks_Morse(options.cmode, kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1)
elif options.wfn.lower() == "ho":
    k_h, k_d = get_ks_HO(kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1)

# Get the highest temperature
T = 1000.0 / WT_HT_list[-1]

if options.wfn.lower() == "ms":
    k_h, k_d = get_ks_Morse(options.cmode, kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1)
elif options.wfn.lower() == "ho":
    k_h, k_d = get_ks_HO(kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1)
"""

# Calculate rates and KIE for various temperatures
Hrate_vsT = []
Drate_vsT = []
KIE_vsT = []

#print('%7s %7s %7s %7s' %('Temper.', 'Rate_H', 'Rate_D', 'KIE'))

for i in xrange(len(T_list1)):

    T = 1000.0 / T_list1[i]

    if options.wfn.lower() == "ms":
        k_h, k_d = get_ks_Morse(options.cmode, kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1)
    elif options.wfn.lower() == "ho":
        k_h, k_d = get_ks_HO(kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1)
    elif options.wfn.lower() == "bs":
        npots = 128
        k_h, k_d = get_ks_bspline(kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, npots, Qm1)

    kie = k_h/k_d

    #print('%7.3e %7.3e %7.3e %7.3e' %(T, k_h, k_d, kie))

    Hrate_vsT.append(k_h)

    if T_list1[i] in T_list2:
        Drate_vsT.append(k_d)
        KIE_vsT.append(kie)

write_file(options.wfn.lower() + '_Hrate_full.txt', T_list1, kh_list, Hrate_vsT)
write_file(options.wfn.lower() + '_Drate.txt', T_list2, kd_list, Drate_vsT)
write_file(options.wfn.lower() + '_KIE.txt', T_list2, KIE_list, KIE_vsT)

quit()

