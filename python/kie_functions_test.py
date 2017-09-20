#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: cal_kie.py
from __future__ import print_function
from math import sqrt
from string_functions import read_list
from numpy import log10
from matplotlib import pyplot as plt
from optparse import OptionParser
from kie_functions import *

###############################################################################
#                                    Tests
###############################################################################

R=2.77
req_ch=1.09
req_oh=0.96

beta_ch=2.0680
beta_oh=2.4420
#for CH, 2900 cm^-1 is an experimental value
#D_ch = (2900.0/(beta_ch*153.467))**2 * ((cmass*hmass)/(cmass+hmass))
#3500 cm^-1 is an experimental value
#D_oh = (3500.0/(beta_oh*153.467))**2 * ((omass*hmass)/(omass+hmass))

D_ch = 77.0
D_oh = 82.0
#beta_ch = (2900.0 / (sqrt(D_ch / ((cmass*hmass)/(cmass+hmass))))) /  153.467
#beta_oh = (3500.0 / (sqrt(D_oh / ((omass*hmass)/(omass+hmass))))) /  153.467
#beta_ch = (2900.0 / (sqrt(D_ch / hmass))) /  153.467
#beta_oh = (3500.0 / (sqrt(D_oh / hmass))) /  153.467

dG0 = -5.4 #kcal/mol
V_el = 4.5 #kcal/mol

coeff = [R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el]
u = 0
v = 0
hmass = 1.0072756064562605

a0 = cal_Suv_trapz(coeff, u, v, hmass)
a1 = cal_Suv_math(coeff, u, v, hmass)
a2 = cal_Suv(coeff, u, v, hmass)
print(a0, a1, a2)

b1 = cal_morse_ene_math(D_ch, hmass, beta_ch, 0)
b2 = cal_morse_ene(D_ch, hmass, beta_ch, 0)
print(b1, b2)

quit()

print('%4s' %'R', end='')
for u in xrange(0, umax+1):
    for v in xrange(0, vmax+1):
        print('%14s' %('S'+str(u)+str(v)), end='')
print('')

for i in range(26, 27):
    R = float(i)*0.1
    print('%4.1f' %R, end='')
    for u in xrange(0, umax+1):
        for v in xrange(0, vmax+1):
                coeff = [R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el]
                Suv = cal_Suv(coeff, u, v, hmass)
                print('%14.6e' %Suv, end='')
    print('')


