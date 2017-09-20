#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: cal_kie.py
from __future__ import print_function
from math import sqrt, gamma, pi
from os import popen, system
from numpy import trapz, inf, real, exp
from scipy.special import eval_genlaguerre as elag
from scipy.special import kn
from scipy.integrate import quad
from scipy.misc import comb

system("module load mathematica/mathematica-11")
system("cp ~/Projs/PCET_code/MathPCETwithET.m .")

###############################################################################
#                                   Wavefunctions
###############################################################################
def cal_wfn(D, m, a, r, re, n):
    #NOT USED!
    """Calculate the wavefunction for Morse potential"""

    #Calculate the constant
    c = (4.184 * 1000.0 / 6.0221409) * 1.660539040
    #                     *=10**-23  *=10**-27  *=10^-10
    c = sqrt(c)/1.0545718
    #           *=10**-34
    c = 0.1 * c
    #           as 10**-35 * 10**34 = 0.1

    # 1 kcal/mol = 4.184*1000.0/6.0221409*10**-23
    # 1 u = 1.660539040*10**-27 #kg
    # hbar = 1.0545718*10**-34  #J*s
    # 1A^-1 = 10**10 m^-1
    # c has the unit of 1

    x = a*r
    xe = a*re
    #m = (m1 * m2) / (m1 + m2)
    lamada = sqrt(2*m*D) * c / a
    z = 2.0 * lamada * exp(xe-x)
    b = lamada-n-0.5
    Nn = sqrt((gamma(n+1.0)*2.0*b)/(gamma(2*lamada-n)))
    psin = Nn * (z**b) * exp(-0.5 * z) * elag(n, 2.0*b, z)
    return psin

def cal_wfn0(D, m, a, r, re):
    #NOT USED!
    """Calculate the ground-state wavefunction for Morse potential"""

    #Calculate the constant
    c = (4.184 * 1000.0 / 6.0221409) * 1.660539040
    #                     *=10**-23  *=10**-27  *=10^-10
    c = sqrt(c)/1.0545718
    #           *=10**-34
    c = 0.1 * c
    #           as 10**-35 * 10**34 = 0.1

    # 1 kcal/mol = 4.184*1000.0/6.0221409*10**-23
    # 1 u = 1.660539040*10**-27 #kg
    # hbar = 1.0545718*10**-34  #J*s
    # 1A^-1 = 10**10 m^-1
    # c has the unit of 1

    x = a*r
    xe = a*re
    #m = (m1 * m2) / (m1 + m2)

    lamada = sqrt(2*m*D) * c / a
    z = 2.0 * lamada * exp(xe-x)
    Nn = sqrt((2.0*lamada-1)/gamma(2.0*lamada))
    #factor = 27.0 / ((2.0 * lamada)**(lamada-0.5) * exp(0.0-lamada) * Nn)
    factor = 1.0
    psi0 = factor * Nn * (z ** (lamada-0.5)) * exp(-0.5 * z)
    return psi0

###############################################################################
#           Overlap integral between two wavefunctions
###############################################################################

def cal_Suv_trapz(coeff, u, v, xmass):
    #NOT USED!
    """Calculate the overlap integral square between donor-H and acceptor-H
       wavefunctions"""

    #global omass, cmass
    R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el = coeff

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

def cal_Suv_Sym(D, m, a, R, u, v):
    #NOT USED!
    """Referred from equation 3.5 in the following paper:
    Anharmonicity effects in atom group transfer processes in condensed phases.
    Journal of the Chemical Society, Faraday Transactions 2:
    Molecular and Chemical Physics 1978, 74, 1690-1701.

    However, there is about 1000.0 times of bug from the numerical results,
    which may comes from a coefficient calculation.
    Also, it assumes we have the same beta values for DH and AH
    """

    #Calculate the constant
    c = (4.184 * 1000.0 / 6.0221409) * 1.660539040 
    #                     *=10^-23  *=10^-27  [*=10^-50 J*kg]
    c = sqrt(c)/1.0545718
    #           *=10^34 sqrt(J*kg)*m/(J*s)
    c = 0.1 * c
    #           as 10^-25 * 10^34 * 10^-10 = 0.1

    lamada = sqrt(2.0*m*D) * c / a #unitless
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

def cal_Suv_math(coeff, u, v, xmass):
    """Calculate the overlap integral square between donor-H and acceptor-H
       wavefunctions based on the Mathmatica code"""

    coeff = [str(i) for i in coeff]
    u = str(u)
    v = str(v)

    #global omass, cmass
    R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el = coeff

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

def cal_Suv(coeff, v1, v2, xmass):
    """Calculate the overlap integral square between donor-H and acceptor-H
       wavefunctions based on python code,
       Corresponding to MorseOverlap in MathPCETwithET.m,
       From equation 20 and Appendix A of International journal of
       quantum chemistry 2002, 88 (2), 280-295."""

    R, D1, beta1, req_ch, D2, beta2, req_oh, dG0, V_el = coeff
    R = R - req_ch - req_oh

    #Calculate the constant
    c1 = (4.184 * 1000.0 / 6.0221409) * 1.660539040
    #                     *=10^-23  *=10^-27  [*=10^-50 J*kg]
    c1 = sqrt(c1)/1.0545718
    #           *=10^34 sqrt(J*kg)*m/(J*s)
    c1 = 0.1 * c1
    #           as 10^-25 * 10^34 * 10^-10 = 0.1

    # 1 kcal/mol = 4.184*1000.0/6.0221409*10**-23
    # 1 u = 1.660539040*10**-27 #kg
    # hbar = 1.0545718*10**-34  #J*s
    # 1A^-1 = 10**10 m^-1
    # c has the unit of 1

    #
    # In front of the summation
    #
    lamada1 = sqrt(2.0*xmass*D1) * c1 / beta1
    lamada2 = sqrt(2.0*xmass*D2) * c1 / beta2
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

def cal_Suv_HO(coeff, v1, v2, freq1, freq2, d, xmass):
    """Corresponding to HarmonicOverlap in MathPCETwithET.m
       Reference: J.-L. Chang, J. Mol. Spectrosc. 232 (2005) 102-104"""
    pass

def cal_Suv_HO_math(coeff, v1, v2, freq1, freq2, d, xmass):
    pass

###############################################################################
#                                   Energies
###############################################################################

def cal_morse_ene(D, m, a, n):
    """Corresponding to MorseEnergy in MathPCETwithET.m,
       Calculate the eigen state energies for Morse potential based on python"""

    #Calculate the constant
    c1 = (2.0 * 4.184 * 1000.0 / 6.0221409 * 1.0) / 1.660539040
    #                         *=10^-23   *=10^20   *=10^27  *=10^24 in total
    c1 = sqrt(c1) * 1.0545718
    #    *=10^12  *=10**-34 *=10**-22 in total
    c1 = c1 * 0.1

    # 1 kcal/mol = 4.184*1000.0/6.022*10**-23
    # 1A^-1 = 10**10 m^-1
    # 1 u = 1.660539040*10**-27 #kg
    # hbar = 1.0545718*10**-34  #J*s
    # c has the unit of 1

    #m = (m1 * m2) / (m1 + m2)
    lamada = sqrt(2.0*m*D) * c1 / a

    c2 = 1.0545718 * sqrt(4184.0/(6.0221409*1.660539040)) * 10.0
    homega = c2 * sqrt(2.0*D*(a**2)/m) * 6.0221409 / 4184.0 #Transfer it back to kcal/mol

    En = homega*((n+0.5)-(1.0/(2.0*lamada))*(n+0.5)**2)
    return En

def cal_morse_ene_math(D, m, a, n):
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

def cal_HO_ene(freq, n):
    """Corresponding to HOEnergy in MathPCETwithET.m"""
    pass

def cal_HO_ene_math(freq, n):
    pass

def norm_prob(ene_list, kb, T):
    """Normalized the probablity based on Boltzmann distribution"""
    ene0 = min(ene_list)
    ene_list2 = [i - ene0 for i in ene_list]
    prob = [exp(-i/(kb*T)) for i in ene_list2]
    probsum = sum(prob)
    prob = [i/probsum for i in prob]
    return prob

###############################################################################
#                              Rate constants
###############################################################################

def cal_kuv_R(coeff, u, v, xmass, cmode, kb, T):
    """Calculate the kuv for R"""

    R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el = coeff

    # Calculate the vibrational state energies for u and v
    # and calculate the overlap integral Suv
    if cmode == 'python':
        Eu = cal_morse_ene(D_ch, xmass, beta_ch, u)
        Ev = cal_morse_ene(D_oh, xmass, beta_oh, v)
        Suv = cal_Suv(coeff, u, v, xmass)
    elif cmode == 'math':
        Eu = cal_morse_ene_math(D_ch, xmass, beta_ch, u)
        Ev = cal_morse_ene_math(D_oh, xmass, beta_oh, v)
        Suv = cal_Suv_math(coeff, u, v, xmass)

    #Calculate the rate constat for kuv for R
    lamb = 13.4 #kcal/mol, re-oragnization energy
    hbar = 1.0545718 #*10**-34  #J*s
    hbar = hbar * 6.0221409 / 4184.0 #*10^-11 kcal/mol*s
    c = (1.0/hbar) * (V_el**2) * sqrt(pi/(lamb*kb*T)) #*10^11 s^-1
    #c = c * (1.0 / hbar) * (V_el**2) * sqrt(pi/(lamb*kb*T))

    dGuv = ((dG0+lamb+Ev-Eu)**2)/(4.0*lamb*kb*T) #unitless
    kuv_R = c * (Suv**2) * exp(-dGuv) #unit 10^11 s^-1
    kuv_R = kuv_R * (10.0**11) #unit s^-1
    return kuv_R

def cal_kR(coeff, umax, vmax, xmass, cmode, kb, T):
    """Calculate the k for R, up to umax and vmax"""

    #global cmass
    R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, dG0, V_el = coeff

    #Calculate the eigen energy states
    if cmode == 'math':
        Eu_list = [cal_morse_ene_math(D_ch, xmass, beta_ch, u) for u in xrange(0, umax+1)]
    elif cmode == 'python':
        Eu_list = [cal_morse_ene(D_ch, xmass, beta_ch, u) for u in xrange(0, umax+1)]

    #Normalize the probabilities
    Pu = norm_prob(Eu_list, kb, T)
    k_R = 0.0
    for u in xrange(0, umax+1):
        for v in xrange(0, vmax+1):
            kuv_R = cal_kuv_R(coeff, u, v, xmass, cmode, kb, T)
            k_R += Pu[u] * kuv_R
    return k_R

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
        print("%7.3e %7.3e %7.3e" %(list1[i], list2[i], list3[i]), file=w_file)
    w_file.close()


