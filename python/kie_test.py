#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
from __future__ import print_function
from kie_functions import *

def test1():
    """Test the harmonic overlap"""

    v1 = 0
    v2 = 0
    freq1 = 2900.0
    freq2 = 3500.0
    d = 1.0

    hoverlap = cal_HO_Suv(v1, v2, freq1, freq2, d, hmass)
    doverlap = cal_HO_Suv(v1, v2, freq1, freq2, d, dmass)

    #print(doverlap)
    print(hoverlap, doverlap)

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
    """Test the Morse overlap"""

    R = 3.05
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

    #beta1 = cal_Morse_beta(D1, freq1)
    #beta2 = cal_Morse_beta(D2, freq2)

    coeff = [R, D1, beta1, req_ch, D2, beta2, req_oh, dG0, V_el]

    v1 = 2
    v2 = 2
    hoverlap = cal_Suv(coeff, v1, v2, hmass)
    doverlap = cal_Suv(coeff, v1, v2, dmass)
    print(hoverlap, doverlap)

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
 
def test3():
    """Test calculating energy for Morse wavefunctions"""

    D = 77.0
    m = hmass
    a = 2.068
    n = 0

    energy = cal_morse_ene(D, m, a, n)
    print(energy)

    n = 1
    energy = cal_morse_ene(D, m, a, n)
    print(energy)

    n = 0
    m = dmass
    energy = cal_morse_ene(D, m, a, n)
    print(energy)

    n = 1
    energy = cal_morse_ene(D, m, a, n)
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

def test4():
    """Test calculating energy for harmonic oscillator wavefunctions"""

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

def test5():
    """Test calculations about Morse potential"""

    D = 77.0
    freq = 2900.0
    m1 = 2.0 * hmass
    m2 = 2.0 * hmass

    a = cal_Morse_beta(D, freq, m1, m2)
    print(a)

    # Should be around 2.068

    freq = cal_Morse_v(D, a, m1, m2)
    print(freq)

    # Should be around 2900.0

def test6():

    m = (cmass * hmass) / (cmass + hmass)
    v = 3000.0 #experimental vibrational frequency of C-H bond, unit cm^-1
    k = cal_HO_k(m, v)
    print(k/2)

    #Results:
    #c3-hc  330.6, GAFF Version 1.81, May 2017

    m = (cmass * cmass) / (cmass + cmass)
    k = 600.0 #c3-c3  300.9, GAFF Version 1.81, May 2017
    v = cal_HO_v(m, k)
    print(v)

    #Results:
    #Experimental C-C vibrational frequency is ~1200 cm^-1

def test7():
    m = 1.0
    n = 1
    freq = 2900.0

    k = cal_HO_k(m, freq)
    print(k)

    energy = cal_HO_ene2(freq, n)
    print(energy)
     
    energy = cal_HO_ene1(k, m, n)
    print(energy)

def test8():
    """Getting the force constants for former PCET work about SLO-1"""
    k = cal_HO_k(110.0, 400.0)
    print(k)
    k = cal_HO_k(14.0, 353.0)
    print(k)
    k = cal_HO_k(100.0, 132.0)
    print(k)
    k = cal_HO_k(10.0, 368.2)
    print(k)


test1()
#test2()
#test3()
#test4()
#test5()
#test6()
#test7()
#test8()

quit()
