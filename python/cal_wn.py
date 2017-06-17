#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: cal_wn.py
from math import exp, sqrt, gamma, pi

def cal_wn(D, a, m1, m2):
    m = (m1 * m2) / (m1 + m2)
    c = 1.0/(2*pi) * sqrt(2 * 4184.0/(6.022 * 1.660539040)) #* 10**12 s^-1
    c = c / 3.0 # 10^4 m-1
    # The speed of light is 3.0 * 10^8 m/s
    c = c * 100.0 # 10^4/10^2
    # 1 m^-1 = 10**-2 cm^-1
    wn = sqrt(D/m) * a * c
    return wn

D = (wn/(a*153.467))**2 * m

#D = 77.0
#a = 2.068
#m1 = 12.0
#m2 = 1.0
#D = 82.0
#a = 2.442
#m1 = 16.0
#m2 = 1.0

#D=76.2329
#a= 0.8828
#m1 = 12.0
#m2 = 1.0

D=82.3209
a= 1.0538
m1 = 16.0
m2 = 1.0

wn = cal_wn(D, a, m1, m2)

print(wn)

