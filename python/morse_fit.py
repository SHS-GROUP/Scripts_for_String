#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: morse_fit.py
from __future__ import print_function
from numpy import array
from math import exp, sqrt, gamma, pi
from matplotlib import pyplot as plt
from matplotlib import rc
from scipy.optimize import minimize
from optparse import OptionParser

parser = OptionParser("Usage: morse_fit.py -s slice_num --fn fig_name")
parser.add_option("-s", dest="slice", type='int',
                  help="Slice number of the 2D PMF")
parser.add_option("--fn", dest="fig_name", type='string',
                  help="Figure name")
(options, args) = parser.parse_args()

#####################################Functions#################################
def read_2d_free_energy(fname):
    free_ene_2d_list = []
    r_fname = open(fname, 'r')
    for line in r_fname:
        line = line.strip('\n')
        line = line.split()
        free_ene_2d_list.append(line)
    r_fname.close()
    free_ene_2d_list = array(free_ene_2d_list)
    return free_ene_2d_list

def read_free_energy_x(fname):
    free_ene_x_list = []
    r_fname = open(fname, 'r')
    for line in r_fname:
        line = line.strip('\n')
        free_ene_x_list.append(float(line))
    r_fname.close()
    binsize = free_ene_x_list[1] - free_ene_x_list[0]
    return free_ene_x_list, binsize

def cal_Uf(x, coeff):

    #U_ch = D_ch * (exp(-2.0*beta_ch*(x-req_ch))-2.0*exp(-1.0*beta_ch*(x-req_ch)))
    #U_oh = D_oh * (exp(-2.0*beta_oh*(req_oh-x))-2.0*exp(-1.0*beta_oh*(req_oh-x))) + delta
    #offset = 0.0 - D_oh + delta
    #U_ch = U_ch - offset
    #U_oh = U_oh - offset

    R, beta_ch, req_ch, beta_oh, req_oh, V_el, delta, offset = coeff

    #D = (2900.0/(a*153.467))**2 * (12.0/13.0) for CH, 2900 cm^-1 is an experimental value
    D_ch = (2900.0/(beta_ch*153.467))**2 * (12.0/13.0)

    #D = (3500.0/(a*153.467))**2 * (16.0/17.0) for OH, 3500 cm^-1 is an experimental value
    D_oh = (3500.0/(beta_oh*153.467))**2 * (16.0/17.0)

    U_ch = D_ch * ((1.0 - exp(-1.0*beta_ch*(R/2.0+x/2.0-req_ch)))**2)  #rCH = R/2 + x/2
    U_oh = D_oh * ((1.0 - exp(-1.0*beta_oh*(R/2.0-x/2.0-req_oh)))**2) + delta  #rOH = R/2 - x/2
    U_ch = U_ch + offset
    U_oh = U_oh + offset
    Uf = 0.5 * (U_ch + U_oh) - 0.5 * sqrt((U_ch - U_oh)**2 + 4*V_el**2)
    return Uf

def cal_Uf_list(coeff):
    global r_list
    Uf_list = []
    for r in r_list:  
        Uf = cal_Uf(r, coeff)
        Uf_list.append(Uf)
    return Uf_list

def gene_weights(num_pts):

    sigma = 3.0

    if num_pts%2 == 0:
        peak_pt = num_pts/2 - 0.5
    elif num_pts%2 == 1:
        peak_pt = num_pts/2

    peak_pt = float(peak_pt)
    weights = []
    for i in xrange(0, num_pts):
        weight = (1.0 / sqrt(2.0 * (sigma**2) * pi)) * exp(-(float(i)-peak_pt)**2/(2.0*sigma**2))
        weights.append(weight)
    return weights

def cal_diff(coeff):

    global Uf_tlist

    Uf_list = cal_Uf_list(coeff)

    # Have weights
    #weights = [1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0]

    #weights = gene_weights(len(Uf_list))
    weights = [1.0 for i in xrange(0, len(Uf_list))]
    diff = [weights[i] * abs(Uf_tlist[i] - Uf_list[i]) for i in xrange(0, len(Uf_list))]

    # Absolute difference
    diff = [abs(Uf_tlist[i] - Uf_list[i]) for i in xrange(0, len(Uf_list))]
    diffsum = sum(diff)/len(diff)

    # Least-squares
    #diff = [abs(Uf_tlist[i] - Uf_list[i]) for i in xrange(0, len(Uf_list))]
    #diff = [i**2 for i in diff]
    #diffsum = sqrt(sum(diff))

    return diffsum

def cal_wn(D, a, m1, m2):
    m = (m1 * m2) / (m1 + m2)
    c = 1.0/(2*math.pi) * sqrt(2 * 4184.0/(6.022 * 1.660539040)) #* 10**12 s^-1
    c = c / 3.0 # 10^4 m-1
    # The speed of light is 3.0 * 10^8 m/s
    c = c * 100.0 # 10^4/10^2
    # 1 m^-1 = 10**-2 cm^-1
    wn = sqrt(D/m) * a * c
    return wn

def plot_morse(coeff, r_list, Uf_tlist, diff, slice_num, figname):

    global Rt

    R, beta_ch, req_ch, beta_oh, req_oh, V_el, delta, offset = coeff

    #D = (2900.0/(a*153.467))**2 * (12.0/13.0) for CH, 2900 cm^-1 is an experimental value
    D_ch = (2900.0/(beta_ch*153.467))**2 * (12.0/13.0)

    #D = (3500.0/(a*153.467))**2 * (16.0/17.0) for OH, 3500 cm^-1 is an experimental value
    D_oh = (3500.0/(beta_oh*153.467))**2 * (16.0/17.0)

    xl = []
    U_chl = []
    U_ohl = []
    Uf_l = []
    #psi1_l = []
    #psi2_l = []

    #offset = 0.0 - D_oh + delta

    overlap = 0.0
    for i in xrange(-130, 200):
        x = float(i) * 0.01
        #U_ch = D_ch * (exp(-2.0*beta_ch*(x-req_ch))-2.0*exp(-1.0*beta_ch*(x-req_ch)))
        #U_oh = D_oh * (exp(-2.0*beta_oh*(req_oh-x))-2.0*exp(-1.0*beta_oh*(req_oh-x))) + delta
        #U_ch = D_ch * (exp(-2.0*beta_ch*(x-req_ch))-2.0*exp(-1.0*beta_ch*(x-req_ch)))
        #U_oh = D_oh * (exp(-2.0*beta_oh*(req_oh-x))-2.0*exp(-1.0*beta_oh*(req_oh-x))) + delta
        #U_ch = U_ch - offset
        #U_oh = U_oh - offset

        U_ch = D_ch * ((1.0 - exp(-1.0*beta_ch*(R/2.0+x/2.0-req_ch)))**2)  #rCH = R/2 + x/2
        U_oh = D_oh * ((1.0 - exp(-1.0*beta_oh*(R/2.0-x/2.0-req_oh)))**2) + delta  #rOH = R/2 - x/2
        U_ch = U_ch + offset
        U_oh = U_oh + offset
        Uf = 0.5 * (U_ch + U_oh) - 0.5 * sqrt((U_ch - U_oh)**2 + 4*V_el**2)

        # Oxygen: 15.9994u
        # Carbon: 12.0107u
        # Hydrogen: 1.007825u
        # Deuterium: 2.014102u
        #psi1 = cal_wfn0(D_ch, 12.0107, 1.007825, beta_ch, req_ch, x)
        #psi2 = cal_wfn0(D_oh, 15.9994, 1.007825, beta_oh, x, req_oh)
        #overlap += min([psi1, psi2])
        #overlap *= 0.01

        xl.append(x)
        U_chl.append(U_ch)
        U_ohl.append(U_oh)
        Uf_l.append(Uf)
        #psi1_l.append(psi1)
        #psi2_l.append(psi2)

    plt.text(-0.9, 62.0, "Slice=%7.4f, Rt=%7.4f, R=%7.4f, UAE=%7.4f" %(slice_num, Rt, R, diff))
    plt.text(-0.9, 58.0, "D(CH)=%7.4f, Beta(CH)=%7.4f, Req(CH)=%7.4f" %(D_ch, beta_ch, req_ch))
    plt.text(-0.9, 54.0, "D(OH)=%7.4f, Beta(OH)=%7.4f, Req(OH)=%7.4f" %(D_oh, beta_oh, req_oh))
    plt.text(-0.9, 50.0, "Vel=%7.4f, Delta=%7.4f, Offset=%7.4f" %(V_el, delta, offset))

    #print('overlap=', overlap)
    #plt.legend

    plt.axis([-1.3, 2.0, -20.0, 65.0])
    plt.plot(r_list, Uf_tlist, 'ko')
    plt.plot(xl, U_chl, 'r-', xl, U_ohl, 'b-', xl, Uf_l, 'g-')
    plt.xlabel(r'CH-OH ($\AA$)')
    plt.ylabel('kcal/mol')
    #plt.plot(xl, psi1_l, 'r--', xl, psi2_l, 'b--')
    plt.savefig(figname, dpi=900)
    plt.close()

def print_coeff(coeff):
    R, beta_ch, req_ch, beta_oh, req_oh, V_el, delta, offset = coeff

    #D = (2900.0/(a*153.467))**2 * (12.0/13.0) for CH, 2900 cm^-1 is an experimental value
    D_ch = (2900.0/(beta_ch*153.467))**2 * (12.0/13.0)

    #D = (3500.0/(a*153.467))**2 * (16.0/17.0) for OH, 3500 cm^-1 is an experimental value
    D_oh = (3500.0/(beta_oh*153.467))**2 * (16.0/17.0)

    print("    R=%7.4f, D(CH)=%7.4f, Beta(CH)=%7.4f, Req(CH)=%7.4f, D(OH)=%7.4f, Beta(OH)=%7.4f, Req(OH)=%7.4f, Vel=%7.4f, Delta=%7.4f, Offset=%7.4f"
             %(R, D_ch, beta_ch, req_ch, D_oh, beta_oh, req_oh, V_el, delta, offset))

#####################################Main Program##############################
tot_r_list, r_binsize = read_free_energy_x('R3-R2_R1.pdf.free_energy.xaxis_crd')
tot_R_list, R_binsize = read_free_energy_x('R3-R2_R1.pdf.free_energy.yaxis_crd')
tot_free_ene = read_2d_free_energy('R3-R2_R1.pdf.free_energy')

#####################################Inital guess##############################
#D_ch = 77.0 #kcal/mol
#beta_ch = 2.068 #A^-1
#req_ch = 1.09 #1.1396 #A

#D_oh = 82.0 #kcal/mol
#beta_oh = 2.442 #A^-1
#req_oh = 0.96 #A #1.5758

#V_el = 4.5 #kcal/mol
#delta = -5.0 #kcal/mol
#offset =  0.0 #kcal/mol

D_ch = 85.7630
beta_ch = 1.9604
req_ch = 1.0103

D_oh = 54.3920
beta_oh = 3.0000
req_oh = 0.8197

V_el = 13.0742
delta = -10.0000
offset = 16.5955

#for i in xrange(0, len(tot_free_ene)):
for i in xrange(options.slice, options.slice+1):
    Rt = tot_R_list[i-1]
    coeff = [Rt, beta_ch, req_ch, beta_oh, req_oh, V_el, delta, offset]
    tot_Uf_tlist = tot_free_ene[:,i-1]
    r_list = []
    Uf_tlist = []
    for j in xrange(0, len(tot_Uf_tlist)):
        if tot_Uf_tlist[j] != 'inf':
            Uf_tlist.append(float(tot_Uf_tlist[j]))
            r_list.append(tot_r_list[j])

    print(Rt)
    print(r_list)
    print(Uf_tlist)

    print('Initial coefficients are:')
    print_coeff(coeff)
    diff = cal_diff(coeff)
    print('    Energy difference is: %7.4f' %diff)
    plot_morse(coeff, r_list, Uf_tlist, diff, options.slice, 'Slice'+ str(options.slice) + '_para_no.pdf')

    # coeff = [Rt, beta_ch, req_ch, beta_oh, req_oh, V_el, delta, offset]
    cons = ( #{'type': 'eq', 'fun': lambda x: x[0] - Rt},
             {'type': 'eq', 'fun': lambda x: x[1] - beta_ch},
             {'type': 'eq', 'fun': lambda x: x[2] - req_ch},
             {'type': 'eq', 'fun': lambda x: x[3] - beta_oh},
             {'type': 'eq', 'fun': lambda x: x[4] - req_oh},
             {'type': 'eq', 'fun': lambda x: x[5] - V_el},
             {'type': 'eq', 'fun': lambda x: x[6] - delta},
             {'type': 'eq', 'fun': lambda x: x[7] - offset}
           )

    #cons = ()

    #D_ch = 77.0 #kcal/mol
    #beta_ch = 2.068 #A^-1
    #req_ch = 1.09 #1.1396 #A

    #D_oh = 82.0 #kcal/mol
    #beta_oh = 2.442 #A^-1
    #req_oh = 0.96 #A #1.5758

    #V_el = 4.5 #kcal/mol
    #delta = -5.0 #kcal/mol
    #offset =  0.0 #kcal/mol

    bnds = ((Rt-R_binsize/2.0, Rt+R_binsize/2.0),
            (1.5, 2.5),
            (0.5, 1.5),
            (2.0, 3.0),
            (0.5, 1.5),
            (0.0, 20.0),
            (-10.0, 10.0),
            (-10.0, 20.0),
           )

    print(bnds)

    #cons = ({'type': 'eq', 'fun': lambda x: x[0]/((2900.0/(x[1]*153.467))**2 * (12.0/13.0))},
           #D_ch and beta_ch to reproduce 2900 cm^-1
    #        {'type': 'eq', 'fun': lambda x: x[3]/((3500.0/(x[4]*153.467))**2 * (16.0/17.0))},
           #D_oh and beta_oh to reproduce 3500 cm^-1
    #       )

    result = minimize(cal_diff, coeff, method='SLSQP', constraints=cons, bounds=bnds)
    #result = minimize(cal_diff, coeff)

    print('Final coefficients are:')
    print_coeff(result.x)
    diff = cal_diff(result.x)
    print('    Energy difference is: %7.4f' %diff)
    plot_morse(result.x, r_list, Uf_tlist, diff, options.slice, options.fig_name)

quit()

