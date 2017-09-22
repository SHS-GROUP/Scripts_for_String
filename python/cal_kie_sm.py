#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
from __future__ import print_function
from math import sqrt
from string_functions import read_list
from numpy import log10
from matplotlib import pyplot as plt
from optparse import OptionParser
from kie_functions import *

parser = OptionParser("Usage: cal_kie_math.py -f R1_free_ene_file -c code_mode -m Morse_mode")
parser.add_option("-f", dest="r1fef", type='string',
                  help="R1 free energy file")
parser.add_option("-c", dest="cmode", type='string',
                  help="Code mode (math or python)")
parser.add_option("-p", dest="pfile", type='string',
                  help="Morse potential file")
(options, args) = parser.parse_args()

###############################################################################
#                                    Constants
###############################################################################
SM_HT_list = [3.5317, 3.4704, 3.4112, 3.3540, 3.2987, 3.2452, 3.1934, 3.1432, 3.0945]
SM_Hrate_list = [70.0, 58.0, 57.0, 54.0, 58.0, 63.0, 56.0, 61.0, 66.0]

SM_DT_list = [3.5317, 3.4704, 3.4112, 3.3540, 3.2987, 3.2452, 3.1934, 3.1432, 3.0945]
SM_Drate_list = [0.18, 0.20, 0.27, 0.3, 0.33, 0.36, 0.41, 0.38, 0.53]

SM_KIE_list = [SM_Hrate_list[i]/SM_Drate_list[i] for i in xrange(len(SM_Hrate_list))]

temp2fit = 4 # 1000.0 / 3.2987 = 303 K

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
umax = 3
vmax = 3
Qm1 = 1.0

#
# Get the parition function parameter
#
T = 1000.0 / SM_HT_list[temp2fit]
k_h, k_d = get_ks(options.cmode, kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1, 1)
Qm1 = SM_Hrate_list[temp2fit] / k_h

# Get the lowest temperature
T = 1000.0 / SM_HT_list[0]
k_h, k_d = get_ks(options.cmode, kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1)

# Get the highest temperature
T = 1000.0 / SM_HT_list[-1]
k_h, k_d = get_ks(options.cmode, kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1)

# Calculate rates and KIE
Hrate_vsT = []
Drate_vsT = []
KIE_vsT = []

for i in xrange(len(SM_HT_list)):
    T = 1000.0 / SM_HT_list[i]
    k_h, k_d = get_ks(options.cmode, kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1)
    kie = k_h/k_d

    Hrate_vsT.append(k_h)
    Drate_vsT.append(k_d)
    KIE_vsT.append(kie)

write_file('Hrate_full.txt', SM_HT_list, SM_Hrate_list, Hrate_vsT)
write_file('Drate.txt', SM_DT_list, SM_Drate_list, Drate_vsT)
write_file('KIE.txt', SM_HT_list, SM_KIE_list, KIE_vsT)

quit()

