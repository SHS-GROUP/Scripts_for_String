#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
from __future__ import print_function
from math import sqrt
from numpy import array, exp, log10
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from kie_functions import *
from matplotlib import rc

###############################################################################
#                            Fitting the Suv along R
###############################################################################

def exp1(x, a, b):
    return a * exp(-b * (x-1.09-0.96))

def exp2(x, a, b, c):
    return a * exp(-b * (x-1.09-0.96) - 1.0/2.0 * (c * (x-1.09-0.96)**2))

def cal_Suv_R():

    hmass = 1.0072756064562605
    dmass = 2.0135514936645316

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

    R_l = []
    S00_hl = []
    S00_dl = []
    kie_l = []

    for i in xrange(41, 80):
        R = float(i) * 0.05
        R_l.append(float(R))
        coeff = [R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el]
        S00_h = cal_Suv(coeff, 5, 5, hmass)
        S00_hl.append(float(S00_h))
        S00_d = cal_Suv(coeff, 5, 5, dmass)
        S00_dl.append(float(S00_d))
        kie = S00_h/S00_d
        kie_l.append(log10(kie))
        print('%7.3e %7.3e' %(R, log10(kie)))

    plt.plot(R_l, kie_l, 'r-')
    plt.savefig('Rvslog10KIE_S55.pdf')
    quit()

    #R_l = [0.0, 1.0, 2.0]
    #S00_hl = [2.0, 0.40, 0.30]
    #S00_dl = [2.0, 0.30, 0.20]

    R_l = array(R_l)
    S00_hl = array(S00_hl)
    S00_dl = array(S00_dl)
    kie_l = array(kie_l)

    popt1, pcov1 = curve_fit(exp1, R_l, S00_hl, p0=[1.0, 20.0], maxfev=10000)
    popt2, pcov2 = curve_fit(exp2, R_l, S00_hl, p0=[1.0, 20.0, 3.0], maxfev=10000)
    popt3, pcov3 = curve_fit(exp1, R_l, S00_dl, p0=[1.0, 30.0], maxfev=10000)
    popt4, pcov4 = curve_fit(exp2, R_l, S00_dl, p0=[1.0, 30.0, 3.0], maxfev=10000)
    plt.plot(R_l, S00_hl, 'bo', label='Hydro')
    plt.plot(R_l, exp1(R_l, *popt1), 'ro', label='HonlyA')
    plt.plot(R_l, exp2(R_l, *popt2), 'go', label='HAandB')
    plt.plot(R_l, S00_hl, 'b*', label='Deter')
    plt.plot(R_l, exp1(R_l, *popt3), 'r*', label='DonlyA')
    plt.plot(R_l, exp2(R_l, *popt4), 'g*', label='DAandB')
    plt.legend()
    plt.savefig('S00.pdf')
    plt.close()

