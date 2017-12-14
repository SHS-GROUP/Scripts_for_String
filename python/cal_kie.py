#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
from __future__ import print_function
from math import sqrt
from string_functions import read_list
from numpy import log10
from matplotlib import pyplot as plt
from optparse import OptionParser
from kie_functions import *

parser = OptionParser("Usage: cal_kie.py -f R1_free_ene_file -c code_mode -p para_file [-w Wavefunction]")

parser.set_defaults(cmode='python', wfn='Morse')

parser.add_option("-f", dest="r1fef", type='string',
                  help="R1 free energy file")
parser.add_option("-c", dest="cmode", type='string',
                  help="Code mode (math or python)")
parser.add_option("-p", dest="pfile", type='string',
                  help="Morse potential file")
parser.add_option("-w", dest="wfn", type='string',
                  help="Wavefunction")
(options, args) = parser.parse_args()

###############################################################################
#                                    Constants
###############################################################################
#WT_HT_list = [3.5971, 3.5336, 3.4722, 3.4130, 3.3557, 3.3003, 3.2468, 3.1949, 3.1447, 3.0960]
#WT_Hrate_list = [213.35, 228.94, 269.96, 275.10, 327.19, 297.11, 329.06, 301.77, 348.08, 327.51]

WT_HT_list = [3.5971, 3.5336, 3.4722, 3.3557, 3.3003, 3.2468, 3.1949, 3.1447]
WT_Hrate_list = [213.35, 228.94, 269.96, 327.19, 297.11, 329.06, 301.77, 348.08]

WT_DT_list = [3.5971, 3.5336, 3.4722, 3.3557, 3.3003, 3.2468, 3.1949, 3.1447]
WT_Drate_list = [2.161, 2.4767, 2.9864, 4.2919, 3.6822, 5.023, 4.2442, 4.0017]

WT_KIE_list = [WT_Hrate_list[i]/WT_Drate_list[i] for i in xrange(len(WT_Hrate_list))]

temp2fit = 4 # The 5th temperature as 1000/3.3 = 303 K
###############################################################################
#                                  Main program
###############################################################################
#
# Read the W(R) list
#
R_list = read_list(options.r1fef, 1)
dR = R_list[1] - R_list[0]
WR_list = read_list(options.r1fef, 2)
para_list = read_para_file(options.pfile)

# para_list for Morse: R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el
# para_list for HO: R, freq1, req1, freq2, req2, dG0, V_el

umax = 3
vmax = 3
Qm1 = 1.0

# Get the parition function parameter
T = 1000.0 / WT_HT_list[temp2fit]

if options.wfn.lower() == "morse":
    k_h, k_d = get_ks_Morse(options.cmode, kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1, 1)
elif options.wfn.lower() == "ho":
    k_h, k_d = get_ks_HO(kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1, 1)

Qm1 = WT_Hrate_list[temp2fit] / k_h

# Get the lowest temperature
T = 1000.0 / WT_HT_list[0]

if options.wfn.lower() == "morse":
    k_h, k_d = get_ks_Morse(options.cmode, kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1)
elif options.wfn.lower() == "ho":
    k_h, k_d = get_ks_HO(kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1)

# Get the highest temperature
T = 1000.0 / WT_HT_list[-1]

if options.wfn.lower() == "morse":
    k_h, k_d = get_ks_Morse(options.cmode, kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1)
elif options.wfn.lower() == "ho":
    k_h, k_d = get_ks_HO(kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1)

# Calculate rates and KIE for various temperatures
Hrate_vsT = []
Drate_vsT = []
KIE_vsT = []

#print('%7s %7s %7s %7s' %('Temper.', 'Rate_H', 'Rate_D', 'KIE'))
for i in xrange(len(WT_HT_list)):
    T = 1000.0 / WT_HT_list[i]

    if options.wfn.lower() == "morse":
        k_h, k_d = get_ks_Morse(options.cmode, kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1)
    elif options.wfn.lower() == "ho":
        k_h, k_d = get_ks_HO(kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1)

    kie = k_h/k_d

    #print('%7.3e %7.3e %7.3e %7.3e' %(T, k_h, k_d, kie))
    Hrate_vsT.append(k_h)
    Drate_vsT.append(k_d)
    KIE_vsT.append(kie)

write_file('Hrate_full.txt', WT_HT_list, WT_Hrate_list, Hrate_vsT)
write_file('Drate.txt', WT_DT_list, WT_Drate_list, Drate_vsT)
write_file('KIE.txt', WT_HT_list, WT_KIE_list, KIE_vsT)

quit()

