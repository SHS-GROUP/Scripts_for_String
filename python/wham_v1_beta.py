#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: wham_v1-beta.py
from __future__ import print_function
from numpy import average, std, arange
from numpy import array, matrix, exp, log, linspace, histogram, transpose
import matplotlib.pyplot as plt
import time
from string_functions import window_data, get_color_dict, write_list, write_xy_lists, read_dat_file, get_rc_data_1D

KbT = 0.596 #Energy unit, which is a global variable

###############################################################################
#                    The first part is about one dimension
###############################################################################

def plot_string_1D(num_imgs, dim, data_dict, ini_wind, final_wind, react_paths, figname):

    color_dict = get_color_dict(dim)
    x = range(1, num_imgs+1) # X axis for plotting
    for i in range(1, dim+1):
        y1 = [] # Y axis for constraints
        y2 = [] # Y axis for plotting
        yerr = [] # Error bar for plotting
        clr = color_dict[i]
        for j in range(ini_wind, final_wind+1):
            y1.append(data_dict[j].equ_dis[i-1])
            avgi = data_dict[j].avg[i-1]
            y2.append(avgi)
            stdi = data_dict[j].std[i-1]
            yerr.append(stdi)
        plt.plot(x, y1, color=clr, linestyle='--')
        plt.plot(x, y2, color=clr, linestyle='-', linewidth=2.0, label=react_paths[i-1])
        plt.errorbar(x, y2, yerr=yerr, ecolor=clr, capsize=2.0, capthick=2.0, barsabove=True)
        plt.legend(bbox_to_anchor=(0.0, 1.02, 1.0, 0.102), loc=3, ncol=2, mode="expand", borderaxespad=0.0)
    plt.ylabel('Distance (Angstrom)')
    plt.xlabel('Number of Images')
    plt.xticks(x)
    plt.savefig(figname)
    plt.close()

def plot_free_ene_1D(num_bins, Prob_RC, xaxis_RC, figname, rc_label):

    global KbT

    # Convert the probability into free energy
    free_ene_RC = [0.0 for i in xrange(num_bins)]
    for i in xrange(num_bins):
        free_ene_RC[i] = -KbT*log(Prob_RC[i])

    min_free_ene = min(free_ene_RC)
    free_ene_RC = [free_ene - min_free_ene for free_ene  in free_ene_RC]

    # Print out the free energy
    write_xy_lists(figname + '.free_energy', xaxis_RC, free_ene_RC)

    # Plot out the free energy
    plt.plot(xaxis_RC, free_ene_RC, 'k-')
    plt.ylabel('Free Energy (kcal/mol)')
    plt.xlabel(rc_label)
    plt.savefig(figname, dpi=900)
    plt.close()

def gene_free_ene_1D(data_dict, num_sims, dim, Fx_old, rc_diml, coefl, num_bins, figname, rc_label):

    global KbT

    data_RC, bins_RC, xaxis_RC = get_rc_data_1D(data_dict, num_sims, rc_diml, coefl, num_bins)

    # Generate the free energy for each bin
    Prob_RC = [0.0 for i in xrange(num_bins)]

    count = 0
    for i in xrange(0, num_sims):
        data_per_sim = len(data_dict[i+1].data)
        for j in xrange(0, data_per_sim):
            each_data_RC = data_RC[count]
            each_count_RC, bins_RCp = histogram(each_data_RC, bins=bins_RC)
            count_indx = list(each_count_RC).index(1)

            #Unbias the free energy
            Ubias = [0.0 for zero_val in xrange(num_sims)]
            for k in xrange(num_sims):
                for ii in xrange(0, dim):
                    equ_dis = data_dict[k+1].equ_dis[ii]
                    constr = data_dict[k+1].constr[ii]
                    samp_dis = data_dict[i+1].data[j,ii]
                    Ubias[k] = Ubias[k] + 0.5 * constr * (samp_dis - equ_dis)**2

            denom = 0.0
            for l in xrange(num_sims):
                data_per_sim2 = len(data_dict[l+1].data)
                denom = denom + data_per_sim2 * exp((Fx_old[l]-Ubias[l])/KbT)

            Prob_RC[count_indx] = Prob_RC[count_indx] + 1.0/denom
            count = count + 1

    plot_free_ene_1D(num_bins, Prob_RC, xaxis_RC, figname, rc_label)

###############################################################################
#                    The second part is about two dimensional
###############################################################################

def get_string_2D(data_dict, ini_image, fin_image, rc_diml1, coefl1, rc_diml2, coefl2):
    crd = []
    crd1 = []
    crd2 = []
    for i in xrange(ini_image, fin_image+1):
        data_val1 = 0.0
        for k in xrange(len(rc_diml1)):
            tmp_data = data_dict[i].equ_dis[rc_diml1[k]-1]
            tmp_data = tmp_data * float(coefl1[k])
            data_val1 = data_val1 + tmp_data
        crd1.append(data_val1)

        data_val2 = 0.0
        for k in xrange(len(rc_diml2)):
            tmp_data = data_dict[i].equ_dis[rc_diml2[k]-1]
            tmp_data = tmp_data * float(coefl2[k])
            data_val2 = data_val2 + tmp_data
        crd2.append(data_val2)

    crd.append(crd1)
    crd.append(crd2)
    return crd

def plot_free_ene_2D(num_bins1, num_bins2, xaxis_RC1, xaxis_RC2, Prob_RC12, figname, xlabel1, xlabel2, colmap, ini_crd=[], fin_crd=[]):

    global KbT

    # Convert the probability into free energy
    free_ene_RC12 = [[0.0 for i in xrange(num_bins2)] for j in xrange(num_bins1)]
    for i in xrange(num_bins1):
        for j in xrange(num_bins2):
            free_ene_RC12[i][j] = -KbT*log(Prob_RC12[i][j])

    max_free_ene = -999.0
    min_free_ene = 999.0
    for i in xrange(num_bins1):
        for j in xrange(num_bins2):
            if free_ene_RC12[i][j] < min_free_ene:
                min_free_ene = free_ene_RC12[i][j]
            if free_ene_RC12[i][j] > max_free_ene:
                max_free_ene = free_ene_RC12[i][j]

    # Print out the free energy
    w_free_ene = open(figname+'.free_energy', 'w')
    for i in xrange(num_bins1):
        for j in xrange(num_bins2):
            free_ene_RC12[i][j] = round(free_ene_RC12[i][j] - min_free_ene, 3)
            print("%-10.3f" %free_ene_RC12[i][j], end=' ', file=w_free_ene)
        print('', file=w_free_ene)
    w_free_ene.close()

    # Print the X and Y axes coordinates
    w_xaxis = open(figname + '.free_energy.xaxis_crd', 'w')
    for i in xrange(num_bins1):
        print("%-10.3f" %xaxis_RC1[i], file=w_xaxis)
    w_xaxis.close()

    # Print the X and Y axes coordinates
    w_yaxis = open(figname + '.free_energy.yaxis_crd', 'w')
    for i in xrange(num_bins2):
        print("%-10.3f" %xaxis_RC2[i], file=w_yaxis)
    w_yaxis.close()

    free_ene_RC12 = transpose(free_ene_RC12)

    # Plot out the free energy
    plt.contourf(xaxis_RC1, xaxis_RC2, free_ene_RC12, 50, cmap=colmap)
    plt.colorbar().set_label('Free Energy (kcal/mol)')
    plt.xlabel(xlabel1)
    plt.ylabel(xlabel2)

    if ini_crd != []:
        plt.scatter(ini_crd[0], ini_crd[1], color='r', marker='o', linewidth=2.0)

    if fin_crd != []:
        plt.scatter(fin_crd[0], fin_crd[1], color='g', marker='o', linewidth=2.0)

    plt.savefig(figname, dpi=900)
    plt.close()

def gene_free_ene_2D(data_dict, num_sims, dim, Fx_old, rc_diml1, coefl1, rc_diml2, coefl2, num_bins1, num_bins2, figname, xlabel1, xlabel2, string_seq, colmap):

    global KbT

    data_RC1, bins_RC1, xaxis_RC1 = get_rc_data_1D(data_dict, num_sims, rc_diml1, coefl1, num_bins1)
    data_RC2, bins_RC2, xaxis_RC2 = get_rc_data_1D(data_dict, num_sims, rc_diml2, coefl2, num_bins2)

    # Generate the free energy for each bin
    Prob_RC12 = [[0.0 for i in xrange(num_bins2)] for j in xrange(num_bins1)]

    count = 0
    for i in xrange(0, num_sims):
        data_per_sim = len(data_dict[i+1].data)
        for j in xrange(0, data_per_sim):
            # For the first dimension
            each_data_RC1 = data_RC1[count]
            each_count_RC1, bins_RCp1 = histogram(each_data_RC1, bins=bins_RC1)
            count_indx1 = list(each_count_RC1).index(1)
            # For the second dimension
            each_data_RC2 = data_RC2[count]
            each_count_RC2, bins_RCp2 = histogram(each_data_RC2, bins=bins_RC2)
            count_indx2 = list(each_count_RC2).index(1)

            #Unbias the free energy
            Ubias = [0.0 for zero_val in xrange(num_sims)]
            for k in xrange(num_sims):
                for ii in xrange(0, dim):
                    equ_dis = data_dict[k+1].equ_dis[ii]
                    constr = data_dict[k+1].constr[ii]
                    samp_dis = data_dict[i+1].data[j,ii]
                    Ubias[k] = Ubias[k] + 0.5 * constr * (samp_dis - equ_dis)**2

            denom = 0.0
            for l in xrange(num_sims):
                data_per_sim2 = len(data_dict[l+1].data)
                denom = denom + data_per_sim2 * exp((Fx_old[l]-Ubias[l])/KbT)

            Prob_RC12[count_indx1][count_indx2] = Prob_RC12[count_indx1][count_indx2] + 1.0/denom
            count = count + 1

    ini_crd = get_string_2D(data_dict, string_seq[0], string_seq[1], rc_diml1, coefl1, rc_diml2, coefl2)
    fin_crd = get_string_2D(data_dict, string_seq[2], string_seq[3], rc_diml1, coefl1, rc_diml2, coefl2)
    plot_free_ene_2D(num_bins1, num_bins2, xaxis_RC1, xaxis_RC2, Prob_RC12, figname, xlabel1, xlabel2, colmap, ini_crd, fin_crd)

###############################################################################
                              #The WHAM iteration
###############################################################################

def wham(dim, react_paths, num_imgs, num_cycles, num_bins, wham_conv):

    start_time1 = time.time()

    global KbT

    # The initial string and final string
    first_i = 1
    first_t = num_imgs
    final_i = num_imgs * (num_cycles-1) + 1
    final_t = final_i + num_imgs - 1
    string_seq = [first_i, first_t, final_i, final_t]

    #############################READ THE DATA FILES###########################

    # All of the following list has dimension of iteration * images
    res_crdl = []
    res_forcel = []
    datal = []

    for i in range(1, num_cycles+1):
        window_dir = './' + str(i) + '/'
        res_crdf = window_dir + 'newconstr.' + str(i) + '.dat'
        res_crd = read_dat_file(res_crdf)
        res_forcef = window_dir + 'force.' + str(i) + '.dat'
        res_force = read_dat_file(res_forcef)
        n1=len(res_crd)
        n2=len(res_force)
        if n1 == n2:
            res_crdl = res_crdl + res_crd
            res_forcel = res_forcel + res_force
            for j in range(1, n1+1):
                datf = window_dir + 'distperframe' + str(j) + '.dat'
                data = read_dat_file(datf)
                datal.append(data)
        else:
            raise ValueError('The number of distances and force constants are '
                             'not the same in the step %d!' %i)

    # Check the three lists are consistent
    if len(res_crdl) == len(res_forcel):
        if len(res_crdl) == len(datal):
            num_sims = len(res_crdl)
        else:
            raise ValueError('The numbers of equlibrium distances and windows are not consistent!')
    else:
        raise ValueError('The number of equlibrium distances and force constants are not consistent!')

    # Generate the data_dict, which contains all the data
    data_dict = {}
    for i in range(0, num_sims):
        window_datai = window_data(res_crdl[i], res_forcel[i], datal[i], dim)
        data_dict[i+1] = window_datai

    #Plot the first and last string
    plot_string_1D(num_imgs, dim, data_dict, first_i, first_t, react_paths, 'First_string.pdf')
    plot_string_1D(num_imgs, dim, data_dict, final_i, final_t, react_paths, 'Last_string.pdf')

    #############################THE WHAM CODE#################################
    Ubiasl = []
    for i in xrange(0, num_sims):
        data_per_sim = len(data_dict[i+1].data)
        for j in xrange(0, data_per_sim):
            for k in xrange(0, num_sims):
                Ubias = 0.0
                for ii in xrange(0, dim): #Along each dimension
                    equ_dis = data_dict[k+1].equ_dis[ii]
                    constr = data_dict[k+1].constr[ii]
                    samp_dis = data_dict[i+1].data[j,ii]
                    Ubias = Ubias + 0.5 * constr * (samp_dis - equ_dis)**2
                Ubias = exp(-Ubias/KbT)
                Ubiasl.append(Ubias)

    #WHAM iteration
    Fx_old = [0.0 for i in xrange(num_sims)]
    Fx_prog = []
    iter = 0
    change = 999.0

    while change > wham_conv:
        # This part is based on the equations 12-15 on the paper
        # of Souaille and Roux on Computer Physics Communications
        # 2001, 135, 40-57
        # Another reference:
        # http://membrane.urmc.rochester.edu/sites/default/files/wham/wham_talk.pdf
        expFx_old = [exp(i/KbT) for i in Fx_old]

        Fx = [0.0 for i in xrange(num_sims)] #Initial free energy
        kk=0

        """
        # The code from the WHAM paper

        for k in xrange(0, num_sims):
            ebfk = 0.0
            for i in xrange(0, num_sims):
                data_per_sim = len(data_dict[i+1].data)
                for l in xrange(0, data_per_sim):
                    bottom = 0.0
                    for j in xrange(0, num_sims):
                        bottom = Ubiasl[kk] * expFx_old[j]
                    ebfk = ebfk + Ubiasl[kk]/bottom
        """             

        for i in xrange(0, num_sims): #Sum over big N
            data_per_sim = len(data_dict[i+1].data) #Sum over little n
            for j in xrange(0, data_per_sim): 
                denom = 0.0 # Obtain the denominator
                Ubiasl_k = []
                for k in xrange(num_sims):
                   denom = denom + expFx_old[k]*Ubiasl[kk]
                   Ubiasl_k.append(Ubiasl[kk])
                   kk = kk + 1
                denom = 1.0/denom
                for k in xrange(num_sims):
                    Fx[k] = Fx[k] + Ubiasl_k[k] * denom

        #Get the updated probability
        Fx = [-KbT*log(Fx[i]) for i in xrange(num_sims)] #Transfer the probability into free energy
        Fx0 = Fx[0] #Normalize the Fx values
        Fx = [Fx[i]-Fx0 for i in xrange(num_sims)]
        Fx_old = Fx #Assign the Fx as Fx_old
 
        if iter <= 1:
            Fx_prog = Fx_prog + Fx
        else:
            Fx_prog = Fx_prog + Fx
            Fx_last = Fx_prog[-num_sims:]
            Fx_last2 = Fx_prog[-2*num_sims:-num_sims]
            Fx_diff = [abs(Fx_last[i] - Fx_last2[i]) for i in range(0, num_sims)]
            change = max(Fx_diff)
        iter = iter + 1

    write_list('Fx.dat', Fx, num_sims)
    write_list('Fx_prog.dat', Fx_prog, num_sims)

    cost_time = time.time() - start_time1

    print("%d Iterations were taken!" %(iter))
    print("It costs %f seconds to finish the WHAM cycle!" %cost_time)
    return data_dict, num_sims, string_seq, Fx_old

###############################################################################
                                   #Examples
###############################################################################
"""
#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
from wham_v1_beta import *

start_time0 = time.time() #Count the time

# 1. Setting variables
dim=3 #Dimension of the reaction coordinates
react_paths = ['C-O', 'O-H', 'C-H'] #Reaction path names
num_imgs=18 #Number of images for the string
num_cycles=5 #Number of iteration for the string calculations
num_bins=20 #Number of bins
wham_conv = 0.001 #WHAM converge creteria

# 2. Do the WHAM cycle
data_dict, num_sims, string_seq, Fx_old = wham(dim, react_paths, num_imgs, num_cycles, num_bins, wham_conv)

# 3. Print the final data
# 3.1. For each dimension
for i in xrange(0, dim):
   plot_dim = i + 1
   figname = 'R' + str(i+1) + '.pdf'
   gene_free_ene_1D(data_dict, num_sims, dim, Fx_old, [plot_dim], [1.0], num_bins, figname, react_paths[i])

# 3.2. For linear combination of 1D system: e.g. R2-R3
gene_free_ene_1D(data_dict, num_sims, dim, Fx_old, [2, 3], [1.0, -1.0], num_bins, 'R2-R3.pdf', 'R2-R3')

# 3.3. For normal 2D system: e.g. R2 vs R3
gene_free_ene_2D(data_dict, num_sims, dim, Fx_old, [2], [1.0], [3], [1.0], num_bins, num_bins, 'R2_R3.pdf', 'R2', 'R3', string_seq)

# 3.4. For 2D system with linear combination: e.g. R1 vs R2-R3
gene_free_ene_2D(data_dict, num_sims, dim, Fx_old, [1], [1.0], [2, 3], [1.0, -1.0], num_bins, num_bins, 'R1_R2-R3.pdf', 'R1', 'R2-R3', string_seq)

cost_time = time.time() - start_time0
print("It costs %f seconds to finish the job!" %cost_time)
quit()
"""
