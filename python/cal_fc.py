#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: cal_wn.py
from __future__ import print_function
from math import exp, sqrt, gamma, pi

def cal_fc(m, omega):

    # 1 kcal/mol = 4.184*1000.0/6.022*10**-23
    # 1 u = 1.660539040*10**-27 #kg
    # hbar = 1.0545718*10**-34  #J*s
    # 1A^-1 = 10**10 m^-1

    # m has unit of kg
    # omega has unit of cm-1

    #omega = omega #* 1.0 * 10**2 # Transfer to m^-1
    #omega = omega #* 3.0 * 10**8 # Unit changes to m^-1 * m/s = s^-1
    m = m * 1.660539040 #*10**-27 #kg
    k = m * omega**2 #unit is kg*s^-2(J/m2) * 9.0*10^-7
    k = k * 1.0 / 4184.0 * 6.022   #* 10^23 # 1 J to 1 kcal/mol
                                   # 1 m^2 = 1*10**20 A^2
                                   # 9.0*10^-7 * 10^3 = 9.0*10^-4
    k = k * 9.0 * (10**-4) * (4 * pi**2)
    return k

#For O-H
#omega = 3500.0 #cm^-1
#m = 16.0/17.0 #u

#For C-H
#omega = 2900.0
#m = 12.0/13.0

omega = 132.8
m = 100.0
k = cal_fc(m, omega)
print(k)

def cal_omega(m, k):
    m = m * 1.660539040 #*10**-27 #kg
    k = k / (1.0 / 4184.0 * 6.022)
    k = k / (9.0 * (10**-4) * (4 * pi**2))
    omega = sqrt(k/m)
    return omega

m = 100.0
k = 126.0
omega = cal_omega(m, k)
print(omega)

