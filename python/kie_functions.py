#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
from __future__ import print_function
from math import sqrt, gamma, pi
from os import popen, system
from numpy import trapz, inf, real, exp
from scipy.special import eval_genlaguerre as elag
from scipy.special import kn, hermite
from numpy.polynomial.hermite import hermval
from numpy import polyfit
from scipy.integrate import quad
from scipy.misc import comb
from string_functions import read_list, read_2d_free_energy
from matplotlib import pyplot as plt
from numpy import polyfit, polyval, linspace

#system("module load mathematica/mathematica-11")
#system("cp ~/Projs/PCET_code/MathPCETwithET.m .")

"""
# Constants from MathPCETwithET.m
kb = 3.16683 * (10.0**-6) # Unit: Hatree/K, MathPCETwithET.m
emass = 9.10938356 * (10**-31) # Unit: kg, from https://en.wikipedia.org/wiki/Electron_rest_mass
Dalton = 1.0 / 5.485799090 * (10**4) # Unit: emass, from https://en.wikipedia.org/wiki/Electron_rest_mass
hmass = 1.0072756064562605 # Unit: Dalton, From MathPCETwithET.m
dmass = 2.0135514936645316 # Unit: Dalton, From MathPCETwithET.m
tmass = 3.0160492 # Unit: Dalton, From MathPCETwithET.m
"""
# Set the lamada
#lamb = 17.7

# Constants
amu = 1.660539040 #Unit: 10^-27 kg; atomic mass unit (Dalton), From https://en.wikipedia.org/wiki/Unified_atomic_mass_unit
hbar = 1.0545718 #Unit: 10^-34 J*s; From https://en.wikipedia.org/wiki/Planck_constant
h = hbar * 2.0 * pi
avg_cons = 6.022140857 #Unit: 10^23/mol, From https://en.wikipedia.org/wiki/Avogadro_constant
boltz_cons = 1.38064852 #Unit: 10^-23 J/K, From https://en.wikipedia.org/wiki/Boltzmann_constant
kcal2j = 4184.0 #Unitless, From https://en.wikipedia.org/wiki/Calorie
lspeed = 2.99792458 #Unit: 10^10 cm/s, From https://en.wikipedia.org/wiki/Speed_of_light
bohr2a = 0.52917721067 #Unit: Angstrom, From https://en.wikipedia.org/wiki/Bohr_radius

kb = boltz_cons * avg_cons #J/K 10^-23 * 10^23 = 1.0
kb = kb / kcal2j #Transfer from J/K to kcal/(mol*K)

omass = 15.99491461956 # Unit: Dalton, From https://en.wikipedia.org/wiki/Oxygen-16
cmass = 12.0 # Unit: Dalton, From https://en.wikipedia.org/wiki/Carbon-12
emass = 5.485799090 * (10**-4) # Unit: Dalton
hmass = 1.007276466879 # Proton not hydrogen, Unit: Dalton, From https://en.wikipedia.org/wiki/Proton
# Hydrogen mass: 1.007825 u, From https://en.wikipedia.org/wiki/Hydrogen_atom
# 1.007825 - 5.485799090 * (10**-4) = 1.007276420091 u, agreed with the proton mass
dmass = 2.01410178 - emass # Deuterium ion not neutral, Unit: Dalton, From https://en.wikipedia.org/wiki/Deuterium
tmass = 3.0160492 - emass # Tritium ion not neutral, Unit: Dalton, From https://en.wikipedia.org/wiki/Tritium

###############################################################################
#                         Wavefunctions for Morse potentials
###############################################################################
def cal_wfn(D, m, a, r, re, n):
    #NOT USED!
    """Calculate the wavefunction for Morse potential"""

    # lamada = sqrt(2.0 * m * D) / (a * hbar)
    # m unit: amu  - kg
    # De unit: kcal/mol - J
    # a unit: A^-1 - m^-1
    # lamada unit: sqrt(kg * kg * m^2 * s^-2) / (m^-1 * kg * m^2 * s^-2 * s)
    #              = kg * m * s^-1 / kg * m * s^-1 = 1

    #Calculate the constant
    c = sqrt (amu * (kcal2j / avg_cons)) / (1.0 * hbar)
    #   sqrt (10^-27 * (1.0 /  10^23  )) / (10^10 * 10^-34) = 0.1
    c = 0.1 * c

    x = a*r
    xe = a*re

    lamada = c * sqrt(2.0*m*D) / a
    z = 2.0 * lamada * exp(xe-x)
    b = lamada-n-0.5
    Nn = sqrt((gamma(n+1.0)*2.0*b)/(gamma(2*lamada-n)))
    psin = Nn * (z**b) * exp(-0.5 * z) * elag(n, 2.0*b, z)
    return psin

def cal_wfn0(D, m, a, r, re):
    #NOT USED!
    """Calculate the ground-state wavefunction for Morse potential"""

    # lamada = sqrt(2.0 * m * D) / (a * hbar)
    # m unit: amu  - kg
    # De unit: kcal/mol - J
    # a unit: A^-1 - m^-1
    # lamada unit: sqrt(kg * kg * m^2 * s^-2) / (m^-1 * kg * m^2 * s^-2 * s)
    #              = kg * m * s^-1 / kg * m * s^-1 = 1

    #Calculate the constant
    c = sqrt (amu * (kcal2j / avg_cons)) / (1.0 * hbar)
    #   sqrt (10^-27 * (1.0 /  10^23  )) / (10^10 * 10^-34) = 0.1
    c = 0.1 * c

    x = a*r
    xe = a*re

    lamada = c * sqrt(2*m*D) / a
    z = 2.0 * lamada * exp(xe-x)
    Nn = sqrt((2.0*lamada-1)/gamma(2.0*lamada))

    # Whether to scale:
    # factor = 27.0 / ((2.0 * lamada)**(lamada-0.5) * exp(0.0-lamada) * Nn)

    factor = 1.0
    psi0 = factor * Nn * (z ** (lamada-0.5)) * exp(-0.5 * z)
    return psi0

###############################################################################
#           Overlap integral between two wavefunctions
###############################################################################

##
## Morse wavefunctions
##

def cal_Morse_Suv_trapz(coeff, u, v, xmass):
    #NOT USED!
    """Calculate the overlap integral square between donor-H and acceptor-H
       wavefunctions"""

    #global omass, cmass
    R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el, lamb = coeff

    #The region to be integrated for x which is rCH-rOH
    xmin = -10.0
    xmax = 10.0
    dxval = 0.01
    xlb = int(xmin / dxval)
    xub = int(xmax / dxval)

    overlap_list = []
    for i in xrange(xlb, xub+1):
        x = float(i) * dxval
        #U_ch = D_ch * ((1.0 - exp(-1.0*beta_ch*(R/2.0+x/2.0-req_ch)))**2)  #rCH = R/2 + x/2
        #U_oh = D_oh * ((1.0 - exp(-1.0*beta_oh*(R/2.0-x/2.0-req_oh)))**2) + delta  #rOH = R/2 - x/2
        psich = cal_wfn(D_ch, xmass, beta_ch, R/2.0+x/2.0, req_ch, u)
        psioh = cal_wfn(D_oh, xmass, beta_oh, R/2.0-x/2.0, req_oh, v)
        overlap_list.append(psich*psioh)
    Suv = trapz(overlap_list, dx=dxval)
    #print('%7.1f %d %d %10.4e' %(R, u, v, Suv))
    return Suv

def cal_Morse_Suv_Sym(D, m, a, R, u, v):
    #NOT USED!
    """Referred from equation 3.5 in the following paper:
    Anharmonicity effects in atom group transfer processes in condensed phases.
    Journal of the Chemical Society, Faraday Transactions 2:
    Molecular and Chemical Physics 1978, 74, 1690-1701.

    However, there is about 1000.0 times of bug from the numerical results,
    which may comes from a coefficient calculation.
    Also, it assumes we have the same beta values for DH and AH
    """

    # lamada = sqrt(2.0 * m * D) / (a * hbar)
    # m unit: amu  - kg
    # De unit: kcal/mol - J
    # a unit: A^-1 - m^-1
    # lamada unit: sqrt(kg * kg * m^2 * s^-2) / (m^-1 * kg * m^2 * s^-2 * s)
    #              = kg * m * s^-1 / kg * m * s^-1 = 1

    #Calculate the constant
    c = sqrt (amu * (kcal2j / avg_cons)) / (1.0 * hbar)
    #   sqrt (10^-27 * (1.0 /  10^23  )) / (10^10 * 10^-34) = 0.1
    c = 0.1 * c

    lamada = c * sqrt(2.0*m*D) / a
    p = 2.0*lamada-1.0
    temp_fac0 = ((p-2.0*u)*(p-2.0*v)*((p+1.0)**(2.0*p))) / (gamma(u+1.0)*gamma(p+1.0-u)*gamma(v+1.0)*gamma(p+1.0-v))

    Suv_pre = 0.0
    for l in xrange(0, u+1):
        for k in xrange(0, v+1):
            temp_fac1 = (-1)**(k+l) * comb(u, l) * comb(v, k)
            temp_fac2 = gamma(p+1.0-u)*gamma(p+1.0-v) / (gamma(p+1.0-u-l)*gamma(p+1.0-v-k))
            temp_fac3 = (p+1.0)**(-(l+k))
            temp_fac4 = (p+1.0)*exp(-a*R/2.0)
            temp_fac5 = exp(-a*R*(p-l-k)/2.0)*kn(abs(l-k),temp_fac4)
            temp_fac6 = temp_fac1 * temp_fac2 * temp_fac3 * temp_fac5
            Suv_pre += temp_fac6

    Suv = 4.0 * temp_fac0 * (Suv_pre**2)
    return Suv

def cal_Morse_Suv_math(coeff, u, v, xmass):
    """Calculate the overlap integral square between donor-H and acceptor-H
       wavefunctions based on the Mathmatica code"""

    coeff = [str(i) for i in coeff]
    u = str(u)
    v = str(v)

    #global omass, cmass
    R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el, lamb = coeff

    if xmass < 2.0:
        a = popen('cal_sij_h.m %s %s %s %s %s %s %s %s %s | awk \'{print $2}\''
                  %(D_ch, beta_ch, D_oh, beta_oh, u, v, R, req_ch, req_oh)).read()
    else:
        a = popen('cal_sij_d.m %s %s %s %s %s %s %s %s %s | awk \'{print $2}\''
                  %(D_ch, beta_ch, D_oh, beta_oh, u, v, R, req_ch, req_oh)).read()

    if '*^' in a:
        a = a.split('*^')
        a = float(a[0]) * 10.0**(int(a[1]))

    return float(a)

def integer(x, k1, k2, y1, y2, a):
    return x**(k1+a*k2-1.0)*exp(-1.0/2.0*(y1*x+y2*(x**a)))

def cal_Morse_Suv(coeff, v1, v2, xmass):
    """Calculate the overlap integral square between donor-H and acceptor-H
       wavefunctions based on python code,
       Corresponding to MorseOverlap in MathPCETwithET.m,
       From equation 20 and Appendix A of
       Lopez, Rivera, Smirnov, and Frank,
       International journal of quantum chemistry 2002, 88 (2), 280-295."""

    R, D1, beta1, req_ch, D2, beta2, req_oh, dG0, V_el, lamb = coeff
    R = R - req_ch - req_oh

    # lamada = sqrt(2.0 * m * D) / (beta * hbar)
    # m unit: amu  - kg
    # De unit: kcal/mol - J
    # beta unit: A^-1 - m^-1
    # lamada unit: sqrt(kg * kg * m^2 * s^-2) / (m^-1 * kg * m^2 * s^-2 * s)
    #              = kg * m * s^-1 / kg * m * s^-1 = 1

    #Calculate the constant
    c = sqrt (amu * (kcal2j / avg_cons)) / (1.0 * hbar)
    #   sqrt (10^-27 * (1.0 /  10^23  )) / (10^10 * 10^-34) = 0.1
    c = 0.1 * c

    #
    # In front of the summation
    #
    lamada1 = c * sqrt(2.0*xmass*D1) / beta1
    lamada2 = c * sqrt(2.0*xmass*D2) / beta2
    j1 = lamada1 - 0.5
    j2 = lamada2 - 0.5
    a = beta2/beta1

    # Factor 0
    N1 = sqrt((gamma(v1+1.0)*2.0*(j1-v1))/(gamma(2*(j1+0.5)-v1)))
    N2 = sqrt((gamma(v2+1.0)*2.0*(j2-v2))/(gamma(2*(j2+0.5)-v2)))
    fac0a = N1*N2*sqrt(a)

    y1 = (2*j1+1.0)*exp(-beta1*(R/2.0))
    y2 = (2*j2+1.0)*exp(-beta2*(R/2.0))

    fac0b = (y1**(j1-v1)) * (y2**(j2-v2))
    fac0 = fac0a*fac0b

    #
    # For the summation
    #
    fcsum = 0.0
    for l1 in xrange(0, v1+1):
        for l2 in xrange(0, v2+1):
            fac1 = (-1)**(l1+l2) / (gamma(l1+1)*gamma(l2+1))
            fac2 = comb(2*j1-v1, v1-l1) * comb(2*j2-v2, v2-l2)
            fac3a = ((2*j1+1.0)**l1) * ((2*j2+1.0)**l2)

            fac3b = exp(-beta1*(R/2.0)*l1) * exp(-beta2*(R/2.0)*l2)

            fac3 = fac3a * fac3b
            fcaux = fac0 * fac1 * fac2 * fac3
            #
            # For integration
            #
            k1 = j1-v1+l1
            k2 = j2-v2+l2

            res = quad(integer, 0.0, inf, args=(k1, k2, y1, y2, -a))
            fcsum += fcaux*res[0]

    #fcsum = real(fcsum)
    return fcsum

def cal_Morse_alpha(coeff, v1, v2, xmass, stepsize):
    """Calculate the first derivative for the Morse wavefunctions overlap
       intergral distance dependence"""

    R, D1, beta1, req_ch, D2, beta2, req_oh, dG0, V_el, lamb = coeff
    coeff_m1 = [R - stepsize, D1, beta1, req_ch, D2, beta2, req_oh, dG0, V_el]
    coeff_p1 = [R + stepsize, D1, beta1, req_ch, D2, beta2, req_oh, dG0, V_el]

    Sm1 = cal_Morse_Suv(coeff_m1, v1, v2, xmass)
    S = cal_Morse_Suv(coeff, v1, v2, xmass)
    Sp1 = cal_Morse_Suv(coeff_p1, v1, v2, xmass)

    der1 = (Sp1-Sm1) / (2.0 * stepsize)
    a = - (1.0/S) * der1
  
    return a

def cal_Morse_ab(coeff, v1, v2, xmass, stepsize):
    """Calculate the first and second derivatives for the overlap intergral
       distance dependence"""

    R, D1, beta1, req_ch, D2, beta2, req_oh, dG0, V_el, lamb = coeff
    coeff_m1 = [R - stepsize, D1, beta1, req_ch, D2, beta2, req_oh, dG0, V_el]
    coeff_p1 = [R + stepsize, D1, beta1, req_ch, D2, beta2, req_oh, dG0, V_el]

    Sm1 = cal_Morse_Suv(coeff_m1, v1, v2, xmass)
    S = cal_Morse_Suv(coeff, v1, v2, xmass)
    Sp1 = cal_Morse_Suv(coeff_p1, v1, v2, xmass)

    der1 = (Sp1-Sm1) / (2.0 * stepsize)
    der2 = (Sp1+Sm1-2.0*S) / (stepsize**2)

    a = - (1.0/S) * der1
    b = - 2.0 * (S * der2 - der1**2) / S**2

    return a, b

##
## Harmonic Oscillators
##

def Hermite(n, x):

    cn = []
    for i in xrange(n):
        cn.append(0.0)
    cn.append(1.0)

    Hnx = hermval(x, cn)
    return Hnx

def cal_HO_Suv(v1, v2, freq1, freq2, d, xmass):
    """Corresponding to HarmonicOverlap in MathPCETwithET.m
       Reference: J.-L. Chang, J. Mol. Spectrosc. 232 (2005) 102-104,
       Mainly from equation 20"""

    # a = m * omega / hbar = m * 2pi * v0 / hbar, v0 in unit as cm^-1
    # Unit: kg * cm^-1 [* cm/s] / (J*s) = (kg/s) / (kg*m^2*s^-2*s) = m^-2 = 10^-20 A^-2
    # Magnitude: 10^-27 * 10^10 / 10^-34 = 10^17 m^-2 = 10^-3 A^-2

    # 2.0 * pi was used because of the omega= 2.0 * pi * v0

    freq1 = freq1 * sqrt(hmass / xmass)
    freq2 = freq2 * sqrt(hmass / xmass)

    c = amu * 2.0 * pi * lspeed / hbar * (10.0**-3)

    a1 = c * xmass * freq1 # a1 has unit as A^-2
    a2 = c * xmass * freq2 # a2 has unit as A^-2

    s = (a1*a2*(d**2))/(a1+a2) # unitless
    A = 2.0 * sqrt(a1*a2) / (a1+a2) # unitless

    #print(s, A)

    fac0 = A*exp(-s)/((2.0**(v1+v2))*gamma(v1+1)*gamma(v2+1))
    fac0 = sqrt(fac0)

    #print(fac0)

    b1 = -(d*a2*sqrt(a1))/(a1+a2) # unitless 
    b2 = d*a1*sqrt(a2)/(a1+a2) # unitless

    #print(b1, b2)

    fcsum = 0.0
    for l1 in xrange(0, v1+1):
        for l2 in xrange(0, v2+1):
            fac1 = comb(v1, l1) * comb(v2, l2)
            #print(v1-l1, b1, v2-l2, b2)
            fac2 = Hermite(v1-l1, b1) * Hermite(v2-l2, b2)
            fac3 = ((2.0*sqrt(a1))**l1) * ((2.0*sqrt(a2))**l2)
            if (l1+l2)%2 == 1:
                fac4 = 0.0
            else:
                fac4 = 1
                K = (l1+l2)/2
                for i in xrange(1, K+1):
                    fac4 = fac4 * (2*i-1)
                fac4 = float(fac4) / ((a1+a2)**K)
            fcsum = fcsum + fac1 * fac2 * fac3 * fac4
            #print(fac1, fac2, fac3, fac4)

    fcsum = fcsum * fac0
    return fcsum

def cal_HO_alpha(v1, v2, freq1, freq2, d, xmass, stepsize):
    """Calculate the first derivative for the Harmonic Oscillator wavefunctions
       overlap intergral distance dependence"""

    dm1 = d - stepsize
    dp1 = d + stepsize

    Sm1 = cal_HO_Suv(v1, v2, freq1, freq2, dm1, xmass)
    S = cal_HO_Suv(v1, v2, freq1, freq2, d, xmass)
    Sp1 = cal_HO_Suv(v1, v2, freq1, freq2, dp1, xmass)

    der1 = (Sp1-Sm1) / (2.0 * stepsize)
    a = - (1.0/S) * der1

    return a

def cal_HO_ab(v1, v2, freq1, freq2, d, xmass, stepsize):
    """Calculate the first and second derivatives for the Harmonic Oscillator
       wavefunctions overlap intergral distance dependence"""

    dm1 = d - stepsize
    dp1 = d + stepsize

    Sm1 = cal_HO_Suv(v1, v2, freq1, freq2, dm1, xmass)
    S = cal_HO_Suv(v1, v2, freq1, freq2, d, xmass)
    Sp1 = cal_HO_Suv(v1, v2, freq1, freq2, dp1, xmass)

    der1 = (Sp1-Sm1) / (2.0 * stepsize)
    der2 = (Sp1+Sm1-2.0*S) / (stepsize**2)

    a = - (1.0/S) * der1
    b = - 2.0 * (S * der2 - der1**2) / S**2

    return a, b

def cal_OH_Suv_math(v1, v2, freq1, freq2, d, xmass):
    pass

###############################################################################
#                         Energies based on wavefunctions
###############################################################################

#
# For Morse potentials
#

def cal_Morse_ene(D, m, a, n):
    """Corresponding to MorseEnergy in MathPCETwithET.m,
       Calculate the eigen state energies for Morse potential based on python"""

    """
    FORMULA 1:

    Equation 43 in
    The Morse oscillator in position space, momentum space, and phase space
    J. P. Dahl, M. Springborg, J. Chem. Phys. 1988, 88(7), 4535

    # En = h * v0 * [ (n + 0.5) - (1 / (2 * lamada)) * (n + 0.5)^2 ]
    # Herein h * v0 = hbar * omega

    # lamada = sqrt(2.0 * m * D) / (a * hbar)
    # m unit: amu  - kg
    # De unit: kcal/mol - J
    # a unit: A^-1 - m^-1
    # hbar unit: J*s
    # lamada unit: sqrt(kg * kg * m^2 * s^-2) / (m^-1 * kg * m^2 * s^-2 * s)
    #              = kg * m * s^-1 / kg * m * s^-1 = 1, so lamada is unitless

    #Calculate the constant
    c = sqrt (amu * (kcal2j / avg_cons)) / (1.0 * hbar)
    #   sqrt (10^-27 * (1.0 /  10^23  )) / (10^10 * 10^-34) = 0.1
    c = c * 0.1

    # lamada = sqrt (2 * m * D) / (a * hbar)
    lamada = c * sqrt(2.0*m*D) / a # unitless

    # Formula: hbar * omega = hbar * sqrt (2 * D * a^2 / m)
    # Magnitude:              10^-34 * sqrt(   10^-23 * 10^20 / 10^-27) = 10^-22
    # Unit: J*s    * sqrt(   J * m^-2 / kg)
    #                = kg * m^2 * s^-1 * sqrt(kg * m^2 * s^-2 * m^-2 / kg)
    #                = kg * m^2 * s^-1 * s^-1
    #                = J

    # v0 = a/(2.0*pi) * sqrt(2.0*D/m)
    c1 = sqrt( kcal2j / (avg_cons * amu))
    hv0 = c1 * (hbar * 2.0 * pi) * a/(2.0*pi) * sqrt(2.0*D/m) * avg_cons / kcal2j * 10.0 #Transfer J back to kcal/mol

    En = hv0 * ((n+0.5)-(1.0/(2.0*lamada))*(n+0.5)**2)
    """

    # FORMULA 2
    # From Wikipedia: https://en.wikipedia.org/wiki/Morse_potential
    # after the line "The eigenenergies in the initial variables have form:"
    # En = h * v0 * (n + 0.5) - [(h * v0 * (n + 0.5))^2 / (4 * D)]

    # v0 = a/(2.0*pi) * sqrt(2.0*D/m)
    c1 = sqrt( kcal2j / (avg_cons * amu))
    hv0 = c1 * (hbar * 2.0 * pi) * a/(2.0*pi) * sqrt(2.0*D/m) * avg_cons / kcal2j * 10.0 #Transfer J back to kcal/mol
    En = hv0 * (n+0.5)-((hv0*(n+0.5))**2)/(4.0*D)
    return En

def cal_Morse_ene_math(D, m, a, n):
    """Calculate the eigen state energies for Morse potential from the Mathmatica code"""

    D = str(D)
    a = str(a)
    n = str(n)

    if m < 2.0:
        #system('cal_morse_ene_h.m %s %s %s' %(D, a, n))
        En = popen('cal_morse_ene_h.m %s %s %s | awk \'{print $2}\'' %(D, a, n)).read()
    else:
        #system('cal_morse_ene_d.m %s %s %s' %(D, a, n))
        En = popen('cal_morse_ene_d.m %s %s %s | awk \'{print $2}\'' %(D, a, n)).read()

    if '*^' in En:
        En = En.split('*^')
        En = float(En[0]) * 10.0**(int(En[1]))
    else:
        En = float(En)

    return En

#
# For Harmonic Oscillators
#

def cal_HO_ene1(k, m, n):
    """Corresponding to HOEnergy in MathPCETwithET.m"""

    # k is the force constant, in unit as kcal/(mol*A^2)
    # m is the reduced mass, in unit as amu
    # n is the state number

    # omega = sqrt(k/m)
    # omega = 2*pi*v0
    # v0 = sqrt(k/m) / (2*pi)  unit: s^-1

    # hv0 = h * sqrt(k/m) / (2*pi) = hbar * sqrt(k/m)
    # Magnitude: 10^-34 * 10^23       * sqrt(1.0 / (10^23 * 10^-27) * 10^20) = 10.0
    # Unit:      kcal/mol * s *       * sqrt(kg * m^2 * s^-2 * m^-2 / kg) =  kcal/mol

    c = hbar * sqrt(kcal2j/(avg_cons * amu)) * avg_cons/kcal2j * 10.0
    hv0 = sqrt(k / m) * c
    energy = hv0 * (n + 0.5)
    return energy

def cal_HO_ene2(v, n):
    """Corresponding to HOEnergy in MathPCETwithET.m"""

    # v is frequency, in unit of cm^-1
    # n the state number

    # v in unit as cm^-1, omega = 2pi * v
    # hv = h * v = hbar * omega
    # Magnitude: 10^-34 * 10^23 * 10^10 = 0.1

    hv = h * v * lspeed * (avg_cons / kcal2j) * 0.1
    energy = hv * (n + 0.5) # in unit of kcal/mol
    return energy

def cal_HO_ene_math(freq, n):
    pass

#
# Normalize the probability based on the energy levels
#

def norm_prob(ene_list, kb, T):
    """Normalized the probablity based on Boltzmann distribution"""

    # Normalize the energy list
    ene0 = min(ene_list)
    ene_list2 = [i - ene0 for i in ene_list]

    # Transfer the energies to probabilities
    prob = [exp(-i/(kb*T)) for i in ene_list2]
    probsum = sum(prob)
    prob = [i/probsum for i in prob]
    return prob

###############################################################################
#               Rate constants for Harmonic Oscillator wavefunctions
###############################################################################

def cal_kuv_R_HO(coeff, u, v, xmass, kb, T):
    """Calculate the kuv for R"""

    R, freq1, req1, freq2, req2, dG0, V_el, lamb = coeff
    R = R - req1 - req2

    # Calculate the vibrational state energies for u and v
    # and calculate the overlap integral Suv
    Eu = cal_HO_ene2(freq1, u)
    Ev = cal_HO_ene2(freq2, v)
    Suv = cal_HO_Suv(u, v, freq1, freq2, R, xmass)

    #Calculate the rate constat for kuv for R

    # Prefac = (V_el**2 / hbar) * sqrt(pi/(lamb*kb*T))
    # Unit: ((kcal/mol)^2 / (J*s)) * sqrt(1.0 / ((kcal/mol)^2))
    #      = (kcal/mol) / (J*s) = J / (J*s) = s^-1
    # Magnitude: (10^-23 / 10^-34) = 10^11 s^-1

    #lamb = 13.4 #kcal/mol, re-oragnization energy
    c = (kcal2j / avg_cons) * (1.0 / hbar)
    Prefac = c * (V_el**2 / hbar) * sqrt(pi/(lamb*kb*T)) #*10^11 s^-1
    dGuv = ((dG0+lamb+Ev-Eu)**2)/(4.0*lamb*kb*T) #unitless

    kuv_R = Prefac * (Suv**2) * exp(-dGuv) #unit 10^11 s^-1
    kuv_R = kuv_R * (10.0**11) #unit s^-1
    return kuv_R

def cal_kR_HO(coeff, umax, vmax, xmass, kb, T):
    """Calculate the k for R, up to umax and vmax"""

    R, freq1, req1, freq2, req2, dG0, V_el, lamb = coeff

    #Calculate the eigen energy states
    Eu_list = [cal_HO_ene2(freq1, u) for u in xrange(0, umax+1)]

    #Normalize the probabilities
    Pu = norm_prob(Eu_list, kb, T)
    k_R = 0.0
    for u in xrange(0, umax+1):
        for v in xrange(0, vmax+1):
            kuv_R = cal_kuv_R_HO(coeff, u, v, xmass, kb, T)
            k_R += Pu[u] * kuv_R
    return k_R

def get_ks_HO(kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1, print_per=0):
    #
    # Get the parition function parameter
    #
    PR_list = norm_prob(WR_list, kb, T)
    k_h_list = []
    k_d_list = []
    for j in xrange(0, len(R_list)):
        R = R_list[j]
        coeff = para_list[j]
        kR_h = cal_kR_HO(coeff, umax, vmax, hmass, kb, T)
        kR_d = cal_kR_HO(coeff, umax, vmax, dmass, kb, T)
        k_h_list.append(kR_h * PR_list[j] * 1.0/float(len(R_list)) * Qm1)
        k_d_list.append(kR_d * PR_list[j] * 1.0/float(len(R_list)) * Qm1)
    k_h = sum(k_h_list)
    k_d = sum(k_d_list)

    #print(k_h_list)
    #print(k_d_list)

    if print_per != 0:
        print("Percentage of %7.1f K" %T)
        print("R", "H_rate", "D_rate")
        for j in xrange(0, len(R_list)):
            print('%6.3f %5.2f %5.2f' %(R_list[j], 100.0 * k_h_list[j]/k_h, 100.0 * k_d_list[j]/k_d))

    return k_h, k_d

###############################################################################
#                     Rate constants for Morse wavefunctions
###############################################################################

def cal_kuv_R_Morse(coeff, u, v, xmass, cmode, kb, T):
    """Calculate the kuv for R"""

    R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el, lamb = coeff

    # Calculate the vibrational state energies for u and v
    # and calculate the overlap integral Suv
    if cmode == 'python':
        Eu = cal_Morse_ene(D_ch, xmass, beta_ch, u)
        Ev = cal_Morse_ene(D_oh, xmass, beta_oh, v)
        Suv = cal_Morse_Suv(coeff, u, v, xmass)
    elif cmode == 'math':
        Eu = cal_Morse_ene_math(D_ch, xmass, beta_ch, u)
        Ev = cal_Morse_ene_math(D_oh, xmass, beta_oh, v)
        Suv = cal_Morse_Suv_math(coeff, u, v, xmass)

    #Calculate the rate constat for kuv for R

    # Prefac = (V_el**2 / hbar) * sqrt(pi/(lamb*kb*T))
    # Unit: ((kcal/mol)^2 / (J*s)) * sqrt(1.0 / ((kcal/mol)^2))
    #      = (kcal/mol) / (J*s) = J / (J*s) = s^-1
    # Magnitude: (10^-23 / 10^-34) = 10^11 s^-1

    #lamb = 13.4 #kcal/mol, re-oragnization energy

    c = (kcal2j / avg_cons) * (1.0 / hbar)
    Prefac = c * (V_el**2 / hbar) * sqrt(pi/(lamb*kb*T)) #*10^11 s^-1
    dGuv = ((dG0+lamb+Ev-Eu)**2)/(4.0*lamb*kb*T) #unitless

    kuv_R = Prefac * (Suv**2) * exp(-dGuv) #unit 10^11 s^-1
    kuv_R = kuv_R * (10.0**11) #unit s^-1
    return kuv_R

def cal_kR_Morse(coeff, umax, vmax, xmass, cmode, kb, T):
    """Calculate the k for R, up to umax and vmax"""

    #global cmass
    R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el, lamb = coeff

    #Calculate the eigen energy states
    if cmode == 'math':
        Eu_list = [cal_Morse_ene_math(D_ch, xmass, beta_ch, u) for u in xrange(0, umax+1)]
    elif cmode == 'python':
        Eu_list = [cal_Morse_ene(D_ch, xmass, beta_ch, u) for u in xrange(0, umax+1)]

    #Normalize the probabilities
    Pu = norm_prob(Eu_list, kb, T)
    k_R = 0.0
    for u in xrange(0, umax+1):
        for v in xrange(0, vmax+1):
            kuv_R = cal_kuv_R_Morse(coeff, u, v, xmass, cmode, kb, T)
            k_R += Pu[u] * kuv_R
    return k_R

def get_ks_Morse(cmode, kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, Qm1, print_per=0):
    #
    # Get the parition function parameter
    #
    PR_list = norm_prob(WR_list, kb, T)
    k_h_list = []
    k_d_list = []
    for j in xrange(0, len(R_list)):
        R = R_list[j]
        coeff = para_list[j]
        kR_h = cal_kR_Morse(coeff, umax, vmax, hmass, cmode, kb, T)
        kR_d = cal_kR_Morse(coeff, umax, vmax, dmass, cmode, kb, T)
        k_h_list.append(kR_h * PR_list[j] * 1.0/float(len(R_list)) * Qm1)
        k_d_list.append(kR_d * PR_list[j] * 1.0/float(len(R_list)) * Qm1)
    k_h = sum(k_h_list)
    k_d = sum(k_d_list)

    #print(k_h_list)
    #print(k_d_list)

    if print_per != 0:
        print("Percentage of %7.1f K" %T)
        print("R", "H_rate", "D_rate")
        for j in xrange(0, len(R_list)):
            print('%6.3f %5.2f %5.2f' %(R_list[j], 100.0 * k_h_list[j]/k_h, 100.0 * k_d_list[j]/k_d))

    return k_h, k_d

###############################################################################
#                     Rate constants for general potentials
###############################################################################

##
## Fourier Grid Hamiltonian method with B-spline
##

def cal_kR_bspline(R, rdatf, outf, mass, kb, T, coeff, umax, vmax, npots):

    # mass is in the unit of electron mass
    Rt, dG0, lamb = coeff

    # Generate a potential plot for target R (Rt)
    Rpl = read_list(rdatf, 1)
    P1l = read_list(rdatf, 4)
    P2l = read_list(rdatf, 6)

    Rpl1 = [(Rp - (Rt - R) / 2.0) for Rp in Rpl]
    Rpl2 = [(Rp + (Rt - R) / 2.0) for Rp in Rpl]

    minval = min(Rpl1)
    maxval = max(Rpl2)
    Rpl_new = linspace(minval, maxval, 20)

    pcoef1 = polyfit(Rpl1, P1l, 10)
    pcoef2 = polyfit(Rpl2, P2l, 10)

    P1_fit = [polyval(pcoef1, i) for i in Rpl_new]
    P2_fit = [polyval(pcoef2, i) for i in Rpl_new]

    datf = outf + '.dat' # The file used for data
    with open(datf, 'w') as f:
        for i in xrange(0, 20):
            print('%5.2f %13.7f %13.7f' %(Rpl_new[i], P1_fit[i], P2_fit[i]), file=f)

    smax = max([umax, vmax])
    system('/share/apps/bspline/bin/fgh_bspline.bin %s %d %10.5f %d > /dev/null ' %(datf, npots, mass, smax))

    overlapf = outf + '_overlaps.dat'
    overlap = read_2d_free_energy(overlapf)
    
    reac_enef = outf + '_reactant_en.dat'
    Eu_list = read_2d_free_energy(reac_enef)
    Eu_list = Eu_list[0]

    prod_enef = outf + '_product_en.dat'
    Ev_list = read_2d_free_energy(prod_enef)
    Ev_list = Ev_list[0]

    # Clean the file
    system("rm %s_*.dat" %outf)

    #
    # Calculate the k at certain R
    #

    # Normalize the probabilities
    Pu = norm_prob(Eu_list, kb, T)
    k_R = 0.0
    V_el = 4.5

    for u in xrange(0, umax):
        for v in xrange(0, vmax):
            Eu = Eu_list[u]
            Ev = Ev_list[v]
            Suv = overlap[u][v]
            #print(Suv)

            #lamb = 13.4 #kcal/mol, re-oragnization energy
            c = (kcal2j / avg_cons) * (1.0 / hbar)
            Prefac = c * (V_el**2 / hbar) * sqrt(pi/(lamb*kb*T)) #*10^11 s^-1
            dGuv = ((dG0+lamb+Ev-Eu)**2)/(4.0*lamb*kb*T) #unitless
            kuv_R = Prefac * (Suv**2) * exp(-dGuv) #unit 10^11 s^-1
            kuv_R = kuv_R * (10.0**11) #unit s^-1
            k_R += Pu[u] * kuv_R

    return k_R

    """"
    Rp1 = Rpl[P1l.index(min(P1l))]
    Rp2 = Rpl[P2l.index(min(P2l))]

    Rchl = []
    Rohl = []
    for Rp in Rpl:
        Rch = R/2.0 + Rp  # Herein Rp = (Rch - Roh)/2.0; R = Rch + Roh
        Roh = R/2.0 - Rp  # So Rch = R/2.0 + Rp; Roh = R/2.0 - Rp
        Rchl.append(Rch)
        Rohl.append(Roh)

    Rp_min = (min(Rchl) - max(Rohl))/2.0
    Rp_max = (max(Rchl) - min(Rohl))/2.0

    plt.plot(Rpl, P1l, 'ro')
    plt.plot(Rpl, P2l, 'bo')
    plt.plot(Rpl1, P1l, 'r^')
    plt.plot(Rpl2, P2l, 'b^')
    plt.plot(Rpl_new, P1_fit, 'r-')
    plt.plot(Rpl_new, P2_fit, 'b-')
    plt.savefig('test.pdf')
    plt.close()

    #plt.plot(Rohl, P2l, 'b')
    #plt.savefig('right.pdf')
    #plt.close()
    """

def get_ks_bspline(kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, npots, Qm1=1.0, print_per=0):
    #
    # Get the parition function parameter
    #

    global emass

    cdft_file = '/home/pengfeil/Projs/SLO_QM/cdft-results-from-alexander/cdft-pt-profiles-RTS.dat'

    k_h_list = []
    k_d_list = []

    PR_list = norm_prob(WR_list, kb, T)
    for j in xrange(0, len(R_list)):
        R = R_list[j]
        coeff = para_list[j]
        kR_h = cal_kR_bspline(2.6, cdft_file, '_temp_bspline', hmass/emass, kb, T, coeff, umax, vmax, npots)
        kR_d = cal_kR_bspline(2.6, cdft_file, '_temp_bspline', dmass/emass, kb, T, coeff, umax, vmax, npots)
                           # (R, rdatf, Rt, outf, hydrogen, kb, T, dG0, umax, vmax, npots):
        k_h_list.append(kR_h * PR_list[j] * 1.0/float(len(R_list)) * Qm1)
        k_d_list.append(kR_d * PR_list[j] * 1.0/float(len(R_list)) * Qm1)
    k_h = sum(k_h_list)
    k_d = sum(k_d_list)

    #print(k_h_list)
    #print(k_d_list)

    if print_per != 0:
        print("Percentage of %7.1f K" %T)
        print("R", "H_rate", "D_rate")
        for j in xrange(0, len(R_list)):
            print('%6.3f %5.2f %5.2f' %(R_list[j], 100.0 * k_h_list[j]/k_h, 100.0 * k_d_list[j]/k_d))

    return k_h, k_d

###############################################################################
#                              Simple calculations
###############################################################################

#
# For Morse potentials
#

def cal_Morse_v(D, a, m1, m2):

    # Formula: v = a / (2*pi) * sqrt(2.0 * D / m)
    # Which is from: https://en.wikipedia.org/wiki/Morse_potential

    # Unit:        m^-1 * sqrt(kg * m^2 * s^-2 / kg) = s^-1
    # Magnitude:   10^10 * sqrt(10^-23 / 10^-27)
    #              = 10^12 s^-1 / [10^10 cm/s] = 100.0 cm^-1

    m = (m1 * m2) / (m1 + m2)
    c = sqrt((kcal2j/avg_cons) / amu) / lspeed * 100.0
    v = c * a / (2.0*pi) * sqrt(2.0*D/m)

    return v

def cal_Morse_beta(D, v, m1, m2):

    # Taking advantage the cal_Morse_v function
    m = (m1 * m2) / (m1 + m2)
    c = sqrt((kcal2j/avg_cons) / amu) / lspeed * 100.0
    a = v  * (2.0*pi) / sqrt(2.0*D/m) / c

    return a

#
# For Harmonic oscillators
#

def cal_HO_k(m1, m2, v):

    # Formula: k = m * omega^2 = m * (2*pi*v)^2
    # Unit:        kg * cm^-2 = kg * [cm^2/s^2] * [(4*pi^2)] * cm^-2
    #                         = kg / s^2 = J / m^2 = (kcal/mol)/(A^2)
    # Magnitude:   10^-27 * [10^20] * 10^23 / 10^20 = 10^-4

    m = (m1 * m2) / (m1 + m2)
    c = amu * (lspeed**2) / kcal2j * avg_cons * (10**-4)
    k = c * m * (2.0 * pi * v)**2

    return k

def cal_HO_v(m1, m2, k):

    # Formula: omega = sqrt(k/m) = 2.0 * pi * v
    # Use the c obtained in function cal_fc
    # Unit of v: cm^-1

    m = (m1 * m2) / (m1 + m2)
    c = amu * (lspeed**2) / kcal2j * avg_cons * (10**-4)
    v = sqrt(k/(c * m)) / (2.0 * pi)

    return v

#
# Return fitted force constant and equilibrium distance
#

def cal_quad_k(rl, enel):

    p = polyfit(rl, enel, 2)
    k = 2.0 * p[0]
    eqr =  -p[1]/(2.0 * p[0])
    return k, eqr

###############################################################################
#                              Reading and Writing
###############################################################################

def read_para_file(fname):
    para_list = []
    readf = open(fname)
    for rline in readf:
        rline = rline.strip('\n')
        line = rline.split()
        line = [float(i) for i in line]
        para_list.append(line)
    readf.close()
    return para_list

def write_file(fname, list1, list2, list3):
    w_file = open(fname, 'w')
    for i in xrange(len(list1)):
        print("%7.4e %7.4e %7.4e" %(list1[i], list2[i], list3[i]), file=w_file)
    w_file.close()


