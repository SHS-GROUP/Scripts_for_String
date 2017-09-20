#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: cal_kie.py
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
DM_HT_list = [3.5317, 3.4704, 3.4112, 3.3540, 3.2987, 3.2452, 3.1934, 3.1432, 3.0945]
#DM_HT_list = [1000.0/i for i in DM_HT_list]
DM_Hrate_list = [70.0, 58.0, 57.0, 54.0, 58.0, 63.0, 56.0, 61.0, 66.0]

DM_DT_list = [3.5317, 3.4704, 3.4112, 3.3540, 3.2987, 3.2452, 3.1934, 3.1432, 3.0945]
#DM_DT_list = [1000.0/i for i in DM_DT_list]
DM_Drate_list = [0.18, 0.20, 0.27, 0.3, 0.33, 0.36, 0.41, 0.38, 0.53]

DM_KIE_list = [DM_Hrate_list[i]/DM_Drate_list[i] for i in xrange(len(DM_Hrate_list))]

temp2fit = 4

T = 300.0 #K
kb = 1.380649 * 6.0221409 #J/K 10^-23 * 10^23 = 1.0
kb = kb / 4184.0 #kcal/(mol*K)

omass = 15.999
cmass = 12.000
hmass = 1.0072756064562605
dmass = 2.0135514936645316

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

#
# Get the parition function parameter
#
T = 1000.0 / DM_HT_list[temp2fit]
PR_list = norm_prob(WR_list, kb, T)
k_h_list = []
k_d_list = []
for j in xrange(0, len(R_list)):
    R = R_list[j]
    coeff = para_list[j]
    kR_h = cal_kR(coeff, umax, vmax, hmass, options.cmode, kb, T)
    kR_d = cal_kR(coeff, umax, vmax, dmass, options.cmode, kb, T)
    k_h_list.append(kR_h * PR_list[j] * 1.0/float(len(R_list)))
    k_d_list.append(kR_d * PR_list[j] * 1.0/float(len(R_list)))
k_h = sum(k_h_list)
k_d = sum(k_d_list)

#print(k_h_list)
#print(k_d_list)

print("Percentage of %7.1f K" %T)
print("R", "H_rate", "D_rate")
for j in xrange(0, len(R_list)):
    print('%6.3f %5.2f %5.2f' %(R_list[j], 100.0 * k_h_list[j]/k_h, 100.0 * k_d_list[j]/k_d))

Qm1 = DM_Hrate_list[temp2fit] / k_h

# Get the lowest temperature

T = 1000.0 / DM_HT_list[0]
PR_list = norm_prob(WR_list, kb, T)
k_h_list = []
k_d_list = []
for j in xrange(0, len(R_list)):
    R = R_list[j]
    coeff = para_list[j]
    kR_h = cal_kR(coeff, umax, vmax, hmass, options.cmode, kb, T)
    kR_d = cal_kR(coeff, umax, vmax, dmass, options.cmode, kb, T)
    k_h_list.append(kR_h * PR_list[j] * 1.0/float(len(R_list)))
    k_d_list.append(kR_d * PR_list[j] * 1.0/float(len(R_list)))
k_h = sum(k_h_list)
k_d = sum(k_d_list)

#print(k_h_list)
#print(k_d_list)

print("Percentage of %7.1f K" %T)
print("R", "H_rate", "D_rate")
for j in xrange(0, len(R_list)):
    print('%6.3f %5.2f %5.2f' %(R_list[j], 100.0 * k_h_list[j]/k_h, 100.0 * k_d_list[j]/k_d))

# Get the highest temperature

T = 1000.0 / DM_HT_list[-1]
PR_list = norm_prob(WR_list, kb, T)
k_h_list = []
k_d_list = []
for j in xrange(0, len(R_list)):
    R = R_list[j]
    coeff = para_list[j]
    kR_h = cal_kR(coeff, umax, vmax, hmass, options.cmode, kb, T)
    kR_d = cal_kR(coeff, umax, vmax, dmass, options.cmode, kb, T)
    k_h_list.append(kR_h * PR_list[j] * 1.0/float(len(R_list)))
    k_d_list.append(kR_d * PR_list[j] * 1.0/float(len(R_list)))
k_h = sum(k_h_list)
k_d = sum(k_d_list)

#print(k_h_list)
#print(k_d_list)

print("Percentage of %7.1f K" %T)
print("R", "H_rate", "D_rate")
for j in xrange(0, len(R_list)):
    print('%6.3f %5.2f %5.2f' %(R_list[j], 100.0 * k_h_list[j]/k_h, 100.0 * k_d_list[j]/k_d))

#
# Calculate rates and KIE
#
Hrate_vsT = []
Drate_vsT = []
KIE_vsT = []

print('%7s %7s %7s %7s' %('Temper.', 'Rate_H', 'Rate_D', 'KIE'))
for i in xrange(len(DM_HT_list)):
    T = 1000.0 / DM_HT_list[i]
    PR_list = norm_prob(WR_list, kb, T)
    #
    #Calculate the KIE value
    #
    k_h_list = []
    k_d_list = []
    for j in xrange(0, len(R_list)):
        R = R_list[j]
        coeff = para_list[j]

        kR_h = cal_kR(coeff, umax, vmax, hmass, options.cmode, kb, T)
        kR_d = cal_kR(coeff, umax, vmax, dmass, options.cmode, kb, T)
        k_h_list.append(kR_h * PR_list[j] * 1.0/float(len(R_list)) * Qm1 )
        k_d_list.append(kR_d * PR_list[j] * 1.0/float(len(R_list)) * Qm1 )

    k_h = sum(k_h_list)
    k_d = sum(k_d_list)
    kie = k_h/k_d
    print('%7.3e %7.3e %7.3e %7.3e' %(T, k_h, k_d, kie))

    Hrate_vsT.append(k_h)
    Drate_vsT.append(k_d)
    KIE_vsT.append(kie)

write_file('Hrate_full.txt', DM_HT_list, DM_Hrate_list, Hrate_vsT)
write_file('Drate.txt', DM_DT_list, DM_Drate_list, Drate_vsT)
write_file('KIE.txt', DM_HT_list, DM_KIE_list, KIE_vsT)

"""
log10KIE_vsT = [log10(i) for i in KIE_vsT]
DM_log10KIE_list = [log10(i) for i in DM_KIE_list]

plt.plot(DM_HT_list, Hrate_vsT, 'rx', label='Theory', markersize=12.0)
plt.plot(DM_HT_list, DM_Hrate_list, 'ro', label='Experiment', markersize=12.0)
plt.legend(bbox_to_anchor=(0.7, 0.2), loc=2, borderaxespad=0.)
#plt.ylim(0.0, 400.0)
plt.xlabel('Temperature (K)')
plt.ylabel(r'Rate($s^{-1}$)')
plt.savefig('HratevsT.pdf', dpi=300.0)
plt.close()

plt.plot(DM_DT_list, Drate_vsT, 'gx', label='Theory', markersize=12.0)
plt.plot(DM_DT_list, DM_Drate_list, 'go', label='Experiment', markersize=12.0)
plt.legend(bbox_to_anchor=(0.7, 0.2), loc=2, borderaxespad=0.)
#plt.ylim(0.0, 6.0)
plt.xlabel('Temperature (K)')
plt.ylabel(r'Rate($s^{-1}$)')
plt.savefig('DratevsT.pdf', dpi=300.0)
plt.close()

plt.plot(DM_DT_list, KIE_vsT, 'bx', label='Theory', markersize=12.0)
plt.plot(DM_DT_list, DM_KIE_list, 'bo', label='Experiment', markersize=12.0)
plt.legend(bbox_to_anchor=(0.7, 0.9), loc=2, borderaxespad=0.)
plt.xlabel('Temperature (K)')
plt.ylabel('KIE')
plt.savefig('KIEvsT.pdf', dpi=300.0)
plt.close()

plt.plot(DM_DT_list, log10KIE_vsT, 'bx', label='Theory', markersize=12.0)
plt.plot(DM_DT_list, DM_log10KIE_list, 'bo', label='Experiment', markersize=12.0)
plt.legend(bbox_to_anchor=(0.7, 0.9), loc=2, borderaxespad=0.)
#plt.ylim(0.0, 5.0)
plt.xlabel('Temperature (K)')
plt.ylabel('log10(KIE)')
plt.savefig('log10KIEvsT.pdf', dpi=300.0)
plt.close()
"""

quit()

