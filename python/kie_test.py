#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
from __future__ import print_function
from kie_functions import *
from string_functions import read_list

#
#
#  About Harmonic Oscillators
#
#

def test1():
    """Test the harmonic overlap"""

    v1 = 0
    v2 = 0
    freq1 = 2900.0
    freq2 = 3500.0
    d = 1.0 # R = 1.09 + 0.96 + d = 3.05 Angstrom

    hoverlap = cal_HO_Suv(v1, v2, freq1, freq2, d, hmass)
    doverlap = cal_HO_Suv(v1, v2, freq1, freq2, d, dmass)
    kie = hoverlap/doverlap

    #print(doverlap)
    print('HO - Hoverlap, Doverlap, KIE:', hoverlap, doverlap, kie)

    # Results from WebPCET server:
    #         http://webpcet.scs.uiuc.edu:8080/webMathematica/webPCET/harmonic.jsp
    # 0, 0
    # hoverlap = 5.131598 * 10^-11
    # doverlap = 2.832093 * 10^-15
    # 1, 1
    # hoverlap = -2.369636 * 10^-9
    # doverlap = -1.860698 * 10^-13
    # 2, 2
    # hoverlap =  5.234826 * 10^-8
    # doverlap =  5.927439 * 10^-12

def test2():
    """Test calculating energy for HO wavefunctions (freqency->energy)"""

    freq = 2900.0
    n = 0
    energy = cal_HO_ene2(freq, n)
    print(energy)
   
    n = 1
    energy = cal_HO_ene2(freq, n)
    print(energy)

    # Results from HOEnergy function in MathPCETwithET.m:
    # freq = 2900.0, n = 0
    # 4.145771013716817
    # freq = 2900.0, n = 1
    # 12.43731304115045

def test2a():
    """Test calculations about different ways of HO frequency->HO energy"""

    m = 1.0
    n = 1
    freq = 2900.0

    k = cal_HO_k(2.0*m, 2.0*m, freq)
    print(k)

    energy = cal_HO_ene2(freq, n)
    print(energy)
     
    energy = cal_HO_ene1(k, m, n)
    print(energy)

def test3():
    """Test calculations about force constant->HO frequency"""

    v = 3000.0 #experimental vibrational frequency of C-H bond, unit cm^-1
    k = cal_HO_k(cmass, hmass, v)
    print(k/2)

    #Results:
    #c3-hc  330.6, GAFF Version 1.81, May 2017

    k = 600.0 #c3-c3  300.9, GAFF Version 1.81, May 2017
    v = cal_HO_v(cmass, cmass, k)
    print(v)

    #Results:
    #Experimental C-C vibrational frequency is ~1200 cm^-1

def test4():
    """Test calculations about HO frequency->force constant"""

    redm = 14.0  # Reduced mass

    v = 250.0 # Frequency
    k = cal_HO_k(2.0*redm, 2.0*redm, v)
    print(k)

    v = 500.0
    k = cal_HO_k(2.0*redm, 2.0*redm, v)
    print(k)

    redm = 100.0

    v = 130.0
    k = cal_HO_k(2.0*redm, 2.0*redm, v)
    print(k)

    v = 300.0
    k = cal_HO_k(2.0*redm, 2.0*redm, v)
    print(k)

    redm = 10.0

    v = 280.0
    k = cal_HO_k(2.0*redm, 2.0*redm, v)
    print(k)

    v = 425.0
    k = cal_HO_k(2.0*redm, 2.0*redm, v)
    print(k)

def test4a():
    """HO frequency->force constant for former PCET work about SLO-1"""

    k = cal_HO_k(2.0*110.0, 2.0*110.0, 400.0)
    print(k)
    k = cal_HO_k(2.0*14.0, 2.0*14.0, 353.0)
    print(k)
    k = cal_HO_k(2.0*100.0, 2.0*100.0, 132.0)
    print(k)
    k = cal_HO_k(2.0*10.0, 2.0*10.0, 368.2)
    print(k)

#
#
#   About Morse potentails
#
#


def test11():
    """Test parameters -> Morse overlap"""

    R = 3.05
    D1 = 77.0
    freq1 = 2900.0
    beta1 = 2.068
    req_ch = 1.09
    D2 = 82.0
    freq2 = 3500.0
    beta2 = 2.442
    req_oh = 0.96
    dG0 = -5.8
    V_el = 4.5
    lamb = 13.4

    #beta1 = cal_Morse_beta(D1, freq1, 2.0*hmass, 2.0*hmass)
    #beta2 = cal_Morse_beta(D2, freq2, 2.0*hmass, 2.0*hmass)

    coeff = [R, D1, beta1, req_ch, D2, beta2, req_oh, dG0, V_el, lamb]

    v1 = 0
    v2 = 0

    hoverlap = cal_Morse_Suv(coeff, v1, v2, hmass)
    doverlap = cal_Morse_Suv(coeff, v1, v2, dmass)
    kie = (hoverlap/doverlap)**2

    print('Morse - Hoverlap, Doverlap, KIE:', hoverlap, doverlap, kie)

    # Results from webPCET server:
    #         http://webpcet.scs.uiuc.edu:8080/webMathematica/webPCET/morse.jsp
    # 0, 0
    # hoverlap = 3.205280 * 10^-7
    # doverlap = 3.298854 * 10^-10
    # 1, 1
    # hoverlap = 3.881105 * 10^-5
    # doverlap = 6.243133 * 10^-8
    # 2, 2
    # hoverlap = 1.601869 * 10^-3
    # doverlap = 4.584307 * 10^-6

def test12():
    """Test parameters->Morse energy"""

    D = 77.0
    m = hmass
    a = 2.068
    n = 0

    energy = cal_Morse_ene(D, m, a, n)
    print(energy)

    n = 1
    energy = cal_Morse_ene(D, m, a, n)
    print(energy)

    n = 0
    m = dmass
    energy = cal_Morse_ene(D, m, a, n)
    print(energy)

    n = 1
    energy = cal_Morse_ene(D, m, a, n)
    print(energy)

    # Results from MorseEnergy function in MathPCETwithET.m:
    #D = 77.0, m = hmass, a = 2.068, n = 0
    # 3.9183572606029364
    #D = 77.0, m = hmass, a = 2.068, n = 1
    # 11.44811603529992
    #D = 77.0, m = dmass, a = 2.068, n = 0
    # 2.7819764392124497
    #D = 77.0, m = dmass, a = 2.068, n = 1
    # 8.192375243310066

def test13():
    """Test frequency->beta and beta->frequency in Morse"""

    D = 77.0
    freq = 2900.0
    m1 = 2.0 * hmass
    m2 = 2.0 * hmass

    a = cal_Morse_beta(D, freq, m1, m2)
    print(a)

    # Should be around 2.068 A^-1

    freq = cal_Morse_v(D, a, m1, m2)
    print(freq)

    # Should be around 2900.0 cm^-1

def test14():
    """Calculating rate constant under 303 K at 3.00 Angstrom with Morse"""

    R = 3.05
    D1 = 77.0
    freq1 = 2900.0
    beta1 = 2.068
    req_ch = 1.09
    D2 = 82.0
    freq2 = 3500.0
    beta2 = 2.442
    req_oh = 0.96
    #dG0 = -5.4
    dG0 = -5.8
    V_el = 1.0
    lamb = 13.4

    #beta1 = cal_Morse_beta(D1, freq1, 2.0*hmass, 2.0*hmass)
    #beta2 = cal_Morse_beta(D2, freq2, 2.0*hmass, 2.0*hmass)

    coeff = [R, D1, beta1, req_ch, D2, beta2, req_oh, dG0, V_el, lamb]

    umax = 2
    vmax = 2

    T = 303.0

    kR_h = cal_kR_Morse(coeff, umax, vmax, hmass, 'python', kb, T)
    kR_d = cal_kR_Morse(coeff, umax, vmax, dmass, 'python', kb, T)

    print('Morse rate constants (H and D), KIE: %7.5e %7.5e %7.5e' %(kR_h, kR_d, kR_h/kR_d))

    # WebPCET result: http://webpcet.scs.uiuc.edu:8080/webMathematica/webPCET/raterfixed.jsp
    # ~ exp(5.11) = 165.67 for Protium at 303 K
    # ~ exp(-6.1) = 0.0022 for Deuterium at 303 K

def test15():
    """Calculating KIE under 303 K at 3.00 Angstrom with Morse"""

    R = 3.00
    D1 = 77.0
    freq1 = 2900.0
    beta1 = 2.068
    req_ch = 1.09
    D2 = 82.0
    freq2 = 3500.0
    beta2 = 2.442
    req_oh = 0.96
    dG0 = -5.4
    V_el = 4.5
    lamb = 13.4

    #beta1 = cal_Morse_beta(D1, freq1, 2.0*hmass, 2.0*hmass)
    #beta2 = cal_Morse_beta(D2, freq2, 2.0*hmass, 2.0*hmass)

    coeff = [R, D1, beta1, req_ch, D2, beta2, req_oh, dG0, V_el, lamb]

    umax = 2
    vmax = 2
 
    T = 303.0

    kR_h = cal_kR_Morse(coeff, umax, vmax, hmass, 'python', kb, T)
    kR_d = cal_kR_Morse(coeff, umax, vmax, dmass, 'python', kb, T)

    print('Morse KIE:', kR_h/kR_d)

    # WebPCET result: http://webpcet.scs.uiuc.edu:8080/webMathematica/webPCET/raterfixed.jsp
    # 7.024060 * 10^4 at 303 K

#
#
#   About general potential
#
#

def test21():
    """Caculate rate constant using general potential under 303 K at 3.05 Angstrom"""

    T = 303.0

    R = 2.6
    dG0 = -5.8
    V_el = 1.0
    lamb = 13.4
    coeff = [R, dG0, V_el, lamb]

    umax = 2
    vmax = 2
    npots = 512
    mass = hmass/emass

    cdft_file = '/home/pengfeil/Projs/SLO_QM/cdft-results-from-alexander/cdft-pt-profiles-RTS-kcal_proton_profiles.dat'
    #'/home/pengfeil/Projs/SLO_QM/cdft-results-from-alexander/cdft-pt-profiles-RTS.dat'

    kR = cal_kR_bspline(2.6, cdft_file, '_temp_bspline', mass, kb, T, coeff, umax, vmax, npots)

    print('General rate constant: %7.5e' %kR)

def test22():
    """Calculate KIE using general potential under 303 K at 3.05 Angstrom"""

    #fname = 'str25+ex3.wt.R3_R1.txt.fit4'

    T = 303.0

    #R_list = read_list(fname, 1)
    #WR_list = read_list(fname, 2)

    R_list = [2.6]
    WR_list = [0.0]

    dG0_list = [-5.8 for i in xrange(len(R_list))]
    V_el_list = [1.0 for i in xrange(len(R_list))]
    lamb_list = [13.4 for i in xrange(len(R_list))]
    para_list = zip(R_list, dG0_list, V_el_list, lamb_list)

    umax = 2
    vmax = 2
    npots = 512

    #cdft_file = '/home/pengfeil/Projs/SLO_QM/cdft-results-from-alexander/cdft-pt-profiles-RTS-kcal_proton_profiles.dat'
    k_h, k_d = get_ks_bspline(kb, T, R_list, WR_list, para_list, umax, vmax, hmass, dmass, npots)

    #print(k_h_list)
    #print(k_d_list)
    print('General potential Hrate, Drate, KIE: %7.5e %7.5e %7.5e' %(k_h, k_d, k_h/k_d))

# Tests about HO
#test1()  # KIE
#test2()
#test2a()
#test3()
#test4()
#test4a()

# Tests about Morse
#test11()  #KIE
#test12()
#test13()
#test14()  #Rate-Morse
#test15()  #KIE-Morse

# Tests about general potential
#test21()  #Rate-bspline
test22()  #KIE-bspline

quit()
