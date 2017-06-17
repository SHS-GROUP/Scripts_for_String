#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: cal_kie.py
from __future__ import print_function
from math import exp, sqrt, gamma, pi
from numpy import trapz
from scipy.special import eval_genlaguerre as elag
from string_functions import read_list
from matplotlib import pyplot as plt

def cal_wfn(D, m1, m2, a, r, re, n):
    """Calculate the wavefunction for Morse potential"""

    #Calculate the constant
    c = (4.184 * 1000.0 / 6.022) * 1.660539040
    #                     *=10**-23  *=10**-27  *=10^-10
    c = sqrt(c)/1.0545718
    #           *=10**-34
    c = 0.1 * c
    #           as 10**-35 * 10**34 = 0.1

    # 1 kcal/mol = 4.184*1000.0/6.022*10**-23
    # 1 u = 1.660539040*10**-27 #kg
    # hbar = 1.0545718*10**-34  #J*s
    # 1A^-1 = 10**10 m^-1
    # c has the unit of 1

    x = a*r
    xe = a*re
    m = (m1 * m2) / (m1 + m2)
    lamada = sqrt(2*m*D) * c / a
    z = 2.0 * lamada * exp(xe-x)
    b = lamada-n-0.5
    Nn = sqrt((gamma(n+1.0)*2.0*b)/(gamma(2*lamada-n)))
    psin = Nn * (z**b) * exp(-0.5 * z) * elag(n, 2.0*b, z)
    return psin

def cal_wfn0(D, m1, m2, a, r, re):
    """Calculate the ground-state wavefunction for Morse potential"""

    #Calculate the constant
    c = (4.184 * 1000.0 / 6.022) * 1.660539040
    #                     *=10**-23  *=10**-27  *=10^-10
    c = sqrt(c)/1.0545718
    #           *=10**-34
    c = 0.1 * c
    #           as 10**-35 * 10**34 = 0.1

    # 1 kcal/mol = 4.184*1000.0/6.022*10**-23
    # 1 u = 1.660539040*10**-27 #kg
    # hbar = 1.0545718*10**-34  #J*s
    # 1A^-1 = 10**10 m^-1
    # c has the unit of 1

    x = a*r
    xe = a*re
    m = (m1 * m2) / (m1 + m2)

    lamada = sqrt(2*m*D) * c / a
    z = 2.0 * lamada * exp(xe-x)
    Nn = sqrt((2.0*lamada-1)/gamma(2.0*lamada))
    #factor = 27.0 / ((2.0 * lamada)**(lamada-0.5) * exp(0.0-lamada) * Nn)
    factor = 1.0
    psi0 = factor * Nn * (z ** (lamada-0.5)) * exp(-0.5 * z)
    return psi0

def cal_morse_ene(D, m1, m2, a, n):
    """Calculate the eigen state energies for Morse potential"""

    #Calculate the constant
    c = (2.0 * 4.184 * 1000.0 / 6.022 * 1.0) / 1.660539040
    #                         *=10**-23  *=10^20   *=10^27  *=10^24 in total
    c = sqrt(c) * 1.0545718
    #    *=10^12  *=10**-34 *=10**-22 in total

    # 1 kcal/mol = 4.184*1000.0/6.022*10**-23
    # 1A^-1 = 10**10 m^-1
    # 1 u = 1.660539040*10**-27 #kg
    # hbar = 1.0545718*10**-34  #J*s
    # c has the unit of 1

    m = (m1 * m2) / (m1 + m2)
    homega = c * sqrt(2.0*D*(a**2)/m)*6.022*10.0/4184.0 #Transfer it back to kcal/mol
    xe = homega/(4.0*D)
    En = homega*((n+0.5)-xe*(n+0.5)**2)
    return En

def norm_prob(ene_list):
    """Normalized the probablity based on Boltzmann distribution"""
    global T

    ene0 = min(ene_list)
    ene_list2 = [i - ene0 for i in ene_list]
    k = 1.38065 * 6.022 #J/K
    k = k/4184.0 #kcal/(mol*K)
    prob = [exp(-i/(k*T)) for i in ene_list2]
    probsum = sum(prob)
    prob = [i/probsum for i in prob]
    return prob

def cal_Suv(coeff, u, v, xmass):
    """Calculate the overlap integral square between donor-H and acceptor-H wavefunctions"""

    global omass, cmass
    R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el = coeff

    #The region to be integrated for x which is rCH-rOH
    xmin = -2.0
    xmax = 2.0
    dxval = 0.01
    xlb = int(xmin / dxval)
    xub = int(xmax / dxval)

    overlap_list = []
    for i in xrange(xlb, xub+1):
        x = float(i) * dxval
        #U_ch = D_ch * ((1.0 - exp(-1.0*beta_ch*(R/2.0+x/2.0-req_ch)))**2)  #rCH = R/2 + x/2
        #U_oh = D_oh * ((1.0 - exp(-1.0*beta_oh*(R/2.0-x/2.0-req_oh)))**2) + delta  #rOH = R/2 - x/2
        psich = cal_wfn(D_ch, cmass, xmass, beta_ch, R/2.0+x/2.0, req_ch, u)
        psioh = cal_wfn(D_oh, omass, xmass, beta_oh, R/2.0-x/2.0, req_oh, v)
        overlap_list.append(psich*psioh)
    Suv = trapz(overlap_list, dx=dxval)
    #print('%7.1f %d %d %10.4e' %(R, u, v, Suv))
    return Suv

def cal_kuv_R(coeff, u, v, xmass):
    """Calculate the kuv for R"""

    global omass, cmass, T
    R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el = coeff

    #Calculate the vibrational state energies for u and v
    Eu = cal_morse_ene(D_ch, cmass, xmass, beta_ch, u)
    Ev = cal_morse_ene(D_oh, omass, xmass, beta_oh, v)

    #Calculate the overlap integral Suv
    Suv = cal_Suv(coeff, u, v, xmass)

    #Calculate the rate constat for kuv for R
    lamb = 13.4 #kcal/mol, re-oragnization energy
    k = 1.38065 * 6.022 #J/K 10^-23 * 10^23 = 1.0
    k = k / 4184.0 #kcal/(mol*K)
    hbar = 1.0545718 #*10**-34  #J*s
    hbar = hbar * 6.022 / 4184.0 #*10^-11 kcal/mol*s
    c = (1.0/hbar) * (V_el**2) * sqrt(pi/(lamb*k*T)) #*10^11 s^-1
    dGuv = ((dG0+lamb+Ev-Eu)**2)/(4.0*lamb*k*T) #unitless
    kuv_R = c * (1.0 / hbar) * (V_el**2) * sqrt(pi/(lamb*k*T)) * (Suv**2) * exp(-dGuv) #unit 10^11 s^-1
    kuv_R = kuv_R * (10.0**11) #unit s^-1
    return kuv_R

def cal_kR(coeff, umax, vmax, xmass):
    """Calculate the k for R, up to umax and vmax"""

    global cmass
    R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el = coeff

    #Calculate the eigen energy states
    Eu_list = [cal_morse_ene(D_ch, cmass, xmass, beta_ch, u) for u in xrange(0, umax+1)]

    #Normalize the probabilities
    Pu = norm_prob(Eu_list)
    k_R = 0.0
    for u in xrange(0, umax+1):
        for v in xrange(0, vmax+1):
            kuv_R = cal_kuv_R(coeff, u, v, xmass)
            k_R += Pu[u] * kuv_R
    return k_R

#
#Constant
#
T = 300.0 #K
omass = 15.9994
cmass = 12.0107
hmass = 1.007825
dmass = 2.014102

#
#Parameters
#
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

dG0 = -5.4 #kcal/mol
V_el = 4.5 #kcal/mol

#
#Calculation the KIE
#
umax = 2
vmax = 2
R_list = read_list('R1.pdf.free_energy', 1)
dR = R_list[1] - R_list[0]
WR_list = read_list('R1.pdf.free_energy', 2)
PR_list = norm_prob(WR_list)

"""
#
#Calculate the KIE value
#
k_h = 0.0
k_d = 0.0

for i in xrange(0, len(R_list)):
    R = R_list[i]
    coeff = [R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el]
    kR_h = cal_kR(coeff, umax, vmax, hmass)
    kR_d = cal_kR(coeff, umax, vmax, dmass)
    k_h += kR_h * PR_list[i] * dR
    k_d += kR_d * PR_list[i] * dR

kie = k_h/k_d
#print("%10.4e" %kie)
"""

#
#For test the overlap intergral calculation
#
print('     R  u v Suv')
for i in range(26, 27):
    R = float(i)*0.1
    for u in xrange(0, umax+1):
        for v in xrange(0, vmax+1):
                coeff = [R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el]
                Suv = cal_Suv(coeff, u, v, hmass)
                print('%7.1f %d %d %10.4e' %(R, u, v, Suv))

quit()

