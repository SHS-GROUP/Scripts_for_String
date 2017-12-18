#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: wham_v1-beta.py
from __future__ import print_function
from numpy import average, std, arange
from numpy import (array, matrix, exp, log, linspace, histogram,
                   transpose, inf, histogramdd)
import matplotlib.pyplot as plt
import time
from string_functions import (window_data, get_color_dict, write_list,
             write_xy_lists, read_dat_file, get_rc_data_1D, read_list, get_data_all)
from random import seed, randint

KbT = 0.596 #Energy unit, which is a global variable

# Reference paper:
# Souaille, M.; Roux, B., Extension to the weighted histogram analysis method:
# combining umbrella sampling with free energy calculations. Computer physics
# communications 2001, 135 (1), 40-57.

###############################################################################
#                    About getting the biased potentials
###############################################################################

def get_Ubiasl(data_dict, num_sims, dim):

    #To get the e^(-beta*Wk(Ri,l)) list (totally has i*l*k terms)
    Ubiasl = []
    for i in xrange(0, num_sims): #Sum over big N for i
        data_per_sim = len(data_dict[i+1].data)
        for l in xrange(0, data_per_sim): #Sum over small n for l
            for k in xrange(0, num_sims): #Sum over i,l for k
                Ubias = 0.0
                for ii in xrange(0, dim): #Along each dimension
                    equ_dis = data_dict[k+1].equ_dis[ii]
                    constr = data_dict[k+1].constr[ii]
                    samp_dis = data_dict[i+1].data[l,ii]
                    Ubias = Ubias + 0.5 * constr * (samp_dis - equ_dis)**2
                Ubias = exp(-Ubias/KbT)
                Ubiasl.append(Ubias)

    # Write into a file
    #w_file = open('expUbiasl.txt', 'w')
    #for i in xrange(len(Ubiasl)):
    #     print("%16.7e" %Ubiasl[i], file=w_file)
    #w_file.close()

    return Ubiasl

###############################################################################
#           The first part is about each image, USUALLY NOT USED
###############################################################################
def plot_wind_avg(x_list, avg_val, xlabel, ylabel, figname):

    w_file = open(figname + '.free_energy', 'w')
    for i in xrange(len(x_list)):
        print("%7.3f %12.7e" %(x_list[i], avg_val[i]), file=w_file)
    w_file.close()

    plt.plot(x_list, avg_val, linestyle='-')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(figname)
    plt.close()

def gene_wind_avg(data_dict, num_sims, dim, expFx, Ubiasl, avg=1, tot_ene=[], ene1=[]):
    """This part is from eqs 19 and 20 from the WHAM paper"""

    if avg not in [1, 2]:
        raise ValueError('avg can only be 1 or 2!')
    elif avg==1:
        data_list = tot_ene
    elif avg==2:
        #Calculate the exp[beta*(U(ri)-U'(ri))]
        ene2 = [tot_ene[i]-ene1[i] for i in xrange(0, len(tot_ene))]
        ene2max = max(ene2)
        ene2min = min(ene2)
        ene2 = [i-(ene2min+ene2max)/2.0 for i in ene2]
        data_list = [exp(i/KbT) for i in ene2]

    avg_val = [0.0 for i in xrange(num_sims)]

    count = 0
    kk = 0
    for i in xrange(0, num_sims): #Sum over big N for i
        data_per_sim = len(data_dict[i+1].data)
        for l in xrange(0, data_per_sim): #Sum over little n for l

            denom = 0.0 #Obtain the denominator, sum over big N for j
            Ubiasl_j = []
            for j in xrange(num_sims):
                data_per_sim2 = len(data_dict[j+1].data)
                denom = denom + float(data_per_sim2) * Ubiasl[kk] * expFx[j]
                Ubiasl_j.append(Ubiasl[kk])
                kk = kk + 1
            denom = 1.0/denom

            for k in xrange(num_sims):
                avg_val[k] = avg_val[k] + data_list[count] * Ubiasl_j[k] * expFx[k] * denom
            count = count + 1

    if avg==2:
        avg_val = [-KbT*log(i) for i in avg_val]
        # Scale
        min_avg_val = min(avg_val)
        avg_val = [i - min_avg_val for i in avg_val]

    return avg_val

###############################################################################
#                    The second part is about one dimension
###############################################################################
def plot_free_ene_1D(num_bins, Prob_RC, xaxis_RC, figname, rc_label, avg=0, count=0):

    global KbT

    if count == 0: # Not printing out the counts
        if avg in [0, 2]:

            # Convert the probability into free energy
            free_ene_RC = [0.0 for i in xrange(num_bins)]
            for i in xrange(num_bins):
                free_ene_RC[i] = -KbT*log(Prob_RC[i])
            # Scale
            min_free_ene = min(free_ene_RC)
            free_ene_RC = [free_ene - min_free_ene for free_ene  in free_ene_RC]

        elif avg == 1:
            free_ene_RC = Prob_RC

        # Print out the free energy
        write_xy_lists(figname + '.free_energy', xaxis_RC, free_ene_RC)

        # Plot out the free energy
        plt.plot(xaxis_RC, free_ene_RC, 'k-')
        plt.ylabel('Free Energy (kcal/mol)')
        plt.xlabel(rc_label)
        plt.savefig(figname, dpi=900)
        plt.close()

    else: # Only for printing out the counts

        # Check the sample number
        n_sample = sum(Prob_RC)
        print('Axis %s has %10d sampled!' %(rc_label, n_sample))

        # Print out the count
        write_xy_lists(figname + '.txt', xaxis_RC, Prob_RC, 1)

        # Plot out the free energy
        plt.plot(xaxis_RC, Prob_RC, 'k-')
        plt.ylabel('Count')
        plt.xlabel(rc_label)
        plt.savefig(figname, dpi=900)
        plt.close()

def gene_free_ene_1D(data_dict, num_sims, dim, expFx, Ubiasl, rc_diml, coefl, num_bins,
                     figname, rc_label, avg=0, tot_ene=[], ene1=[]):
    """The part is from eqs 10 and 11 in the WHAM paper"""

    num_binsX, data_RC, bins_RC, xaxis_RC = get_rc_data_1D(data_dict, num_sims, rc_diml, coefl, num_bins)

    if avg == 0: # No average to calculate
        pass
    elif avg == 1: # Only average to calculate
        if tot_ene == []:
            raise ValueError('No value list provided to calculate the average value!')
        elif tot_ene != []:
            data_list = tot_ene
            if ene1 != []:
                raise ValueError('Provide two list values to calculate the average '
                                 'value, only one list is needed!')
    elif avg == 2: # About the free energy partitioning
        if tot_ene != [] and ene1 != []:
            if len(tot_ene) == len(ene1):
                #Calculate the exp[beta*(U(ri)-U'(ri))]
                ene2 = [tot_ene[i]-ene1[i] for i in xrange(0, len(tot_ene))]
                ene2min = min(ene2)
                ene2 = [i-ene2min-300.0 for i in ene2]
                data_list = [exp(i/KbT) for i in ene2]
            else:
                raise ValueError('The length of total energy list and paritioning '
                                 'energy list are not the same!')
        elif tot_ene != []:
            raise ValueError('Only total energy list provided, you need to '
                             'provide paritioning energy list as well!')
        elif ene1 != []:
            raise ValueError('Only paritioning energy list provided, you need to '
                             'provide total energy list as well!')

    Prob_RCN = [0 for i in xrange(num_binsX)]
    Prob_RC0 = [0.0 for i in xrange(num_binsX)]
    Prob_RC1 = [0.0 for i in xrange(num_binsX)]

    count = 0
    kk = 0

    if avg==0:
        for i in xrange(0, num_sims): #Sum over big N for i
            data_per_sim = len(data_dict[i+1].data)
            for l in xrange(0, data_per_sim): #Sum over little n for l
 
                #Calculate the probability for each point
                each_data_RC = data_RC[count]
                each_count_RC, bins_RCp = histogram(each_data_RC, bins=bins_RC)
                count_indx = list(each_count_RC).index(1)
 
                #Get the biased potentials
                Ubiasl_j = []
                for j in xrange(num_sims):
                    Ubiasl_j.append(Ubiasl[kk])
                    kk = kk + 1
 
                #Obtain the denominator
                denom = 0.0
                for k in xrange(num_sims):
                    data_per_sim2 = len(data_dict[k+1].data)
                    denom = denom + float(data_per_sim2) * Ubiasl_j[k] * expFx[k]
 
                Prob_RCN[count_indx] = Prob_RCN[count_indx] + 1
                Prob_RC0[count_indx] = Prob_RC0[count_indx] + 1.0/denom
                count = count + 1

        plot_free_ene_1D(num_binsX, Prob_RCN, xaxis_RC, figname + '.count.pdf', rc_label, 0, 1)
        plot_free_ene_1D(num_binsX, Prob_RC0, xaxis_RC, figname, rc_label, avg)

    else:
        for i in xrange(0, num_sims): #Sum over big N for i
            data_per_sim = len(data_dict[i+1].data)
            for l in xrange(0, data_per_sim): #Sum over little n for l
 
                # Calculate the probability for each point
                each_data_RC = data_RC[count]
                each_count_RC, bins_RCp = histogram(each_data_RC, bins=bins_RC)
                count_indx = list(each_count_RC).index(1)
 
                # Get the biased potentials
                Ubiasl_j = []
                for j in xrange(num_sims):
                    Ubiasl_j.append(Ubiasl[kk])
                    kk = kk + 1
 
                # Obtain the denominator
                denom = 0.0
                for k in xrange(num_sims):
                    data_per_sim2 = len(data_dict[k+1].data)
                    denom = denom + float(data_per_sim2) * Ubiasl_j[k] * expFx[k]

                Prob_RCN[count_indx] = Prob_RCN[count_indx] + 1
                Prob_RC0[count_indx] = Prob_RC0[count_indx] + 1.0/denom
                Prob_RC1[count_indx] = Prob_RC1[count_indx] + data_list[count]/denom
                count = count + 1

        Prob_RC2 = [Prob_RC1[i]/Prob_RC0[i] for i in xrange(num_binsX)]
        plot_free_ene_1D(num_binsX, Prob_RCN, xaxis_RC, figname + '.count.pdf', rc_label, 0, 1)
        plot_free_ene_1D(num_binsX, Prob_RC2, xaxis_RC, figname, rc_label, avg)

###############################################################################
#                    The third part is about two dimensional
###############################################################################
def get_string_2D(data_dict, ini_image, fin_image, rc_diml1, coefl1, rc_diml2, coefl2, fname):
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

    w_fname = open(fname, 'w')
    for i in xrange(0, len(crd1)):
        print("%7.4f %7.4f" %(crd1[i], crd2[i]), file=w_fname)
    w_fname.close()

    return crd

def plot_string_free_energy(free_ene_RC12, fin_crd, xaxis_RC1, xaxis_RC2, num_bins1, num_bins2, figname):

    free_ene_RC12 = transpose(free_ene_RC12)

    # Plot out the 1D free energy along the final string
    str_free_ene = []
    for i in xrange(0, len(fin_crd[0])):
        ind_dict = {}
        k = 0
        x = fin_crd[0][i]
        y = fin_crd[1][i]
        crd_disl = []
        for i in xrange(num_bins1):
            for j in xrange(num_bins2):
                crd_dis = (xaxis_RC1[i]-x)**2 + (xaxis_RC2[j]-y)**2
                crd_disl.append(crd_dis)
                ind_dict[k] = (i, j)
                k = k + 1
        k = min(crd_disl)
        k = crd_disl.index(k)
        ind_2d = (ind_dict[k][0], ind_dict[k][1])
        free_ene = free_ene_RC12[ind_2d[0]][ind_2d[1]]
        str_free_ene.append(free_ene)
        #print(min(crd_disl), ind_dict[k][0], ind_dict[k][1], free_ene, xaxis_RC1[ind_dict[k][0]], xaxis_RC2[ind_dict[k][1]])

    w_str_free_ene = open(figname + '.free_ene', 'w')
    for i in xrange(0, len(str_free_ene)):
        print("%-10.3f" %str_free_ene[i], file=w_str_free_ene)
    w_str_free_ene.close()

    plt.plot(range(1, len(str_free_ene)+1), str_free_ene, 'k-', linewidth=2.0)
    plt.xticks(range(1, len(str_free_ene)+1))
    plt.xlabel('Number of Images')
    plt.ylabel('Free Energy (kcal/mol)')
    plt.savefig(figname + '.free_ene.pdf', dpi=900)
    plt.close()

def plot_free_ene_2D(num_bins1, num_bins2, xaxis_RC1, xaxis_RC2, Prob_RC12, figname,
                     xlabel1, xlabel2, colmap, ini_crd=[], fin_crd=[], avg=0, count=0):

    global KbT

    if count == 0: # Not print out the counts
        if avg in [0, 2]:
            # Convert the probability into free energy
            free_ene_RC12 = [[0.0 for i in xrange(num_bins2)] for j in xrange(num_bins1)]
            for i in xrange(num_bins1):
                for j in xrange(num_bins2):
                    free_ene_RC12[i][j] = -KbT*log(Prob_RC12[i][j])

            # Normalize the free energy
            min_free_ene = 99999999999.0
            for i in xrange(num_bins1):
                for j in xrange(num_bins2):
                    if free_ene_RC12[i][j] < min_free_ene:
                        min_free_ene = free_ene_RC12[i][j]
  
            for i in xrange(num_bins1):
                for j in xrange(num_bins2):
                    free_ene_RC12[i][j] = round(free_ene_RC12[i][j] - min_free_ene, 3)

        elif avg == 1:
            free_ene_RC12 = Prob_RC12

        ############################
        # Print out
        ############################
        # Print out the free energy
        w_free_ene = open(figname+'.free_energy', 'w')
        for i in xrange(num_bins1):
            for j in xrange(num_bins2):
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
 
        ############################
        # Plotting
        ############################
        # Transpose for plotting
        free_ene_RC12 = transpose(free_ene_RC12)
 
        # Plot out the 2D free energy
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
 
        plot_string_free_energy(free_ene_RC12, ini_crd, xaxis_RC1, xaxis_RC2,
                                num_bins1, num_bins2, figname + '.inistr')
        plot_string_free_energy(free_ene_RC12, fin_crd, xaxis_RC1, xaxis_RC2,
                                num_bins1, num_bins2, figname + '.finstr')

    else: # Only for printing out the counts

        n_sample = 0
        for i in xrange(num_bins1):
            for j in xrange(num_bins2):
                n_sample = n_sample + Prob_RC12[i][j]
        print('Axes %s %s has %10d sampled!' %(xlabel1, xlabel2, n_sample))

        ############################
        # Print out
        ############################
        # Print out the 2D counts
        w_free_ene = open(figname+'.txt', 'w')
        for i in xrange(num_bins1):
            for j in xrange(num_bins2):
                print("%10d" %Prob_RC12[i][j], end=' ', file=w_free_ene)
            print('', file=w_free_ene)
        w_free_ene.close()

        ############################
        # Plotting
        ############################
        # Transpose for plotting
        Prob_RC12 = transpose(Prob_RC12)

        # Plot out map of the 2D counts
        plt.contourf(xaxis_RC1, xaxis_RC2, Prob_RC12, 50, cmap=colmap)
        plt.colorbar().set_label('Count')
        plt.xlabel(xlabel1)
        plt.ylabel(xlabel2)

        if ini_crd != []:
            plt.scatter(ini_crd[0], ini_crd[1], color='r', marker='o', linewidth=2.0)

        if fin_crd != []:
            plt.scatter(fin_crd[0], fin_crd[1], color='g', marker='o', linewidth=2.0)

        plt.savefig(figname, dpi=900)
        plt.close()

        plot_string_free_energy(Prob_RC12, ini_crd, xaxis_RC1, xaxis_RC2,
                                num_bins1, num_bins2, figname + '.inistr')
        plot_string_free_energy(Prob_RC12, fin_crd, xaxis_RC1, xaxis_RC2,
                                num_bins1, num_bins2, figname + '.finstr')

def gene_free_ene_2D(data_dict, num_sims, dim, expFx, Ubiasl, rc_diml1, coefl1,
                     rc_diml2, coefl2, num_bins1, num_bins2, figname, xlabel1, xlabel2,
                     string_seq, colmap, avg=0, tot_ene=[], ene1=[]):

    num_binsX, data_RC1, bins_RC1, xaxis_RC1 = get_rc_data_1D(data_dict, num_sims, rc_diml1, coefl1, num_bins1)
    num_binsY, data_RC2, bins_RC2, xaxis_RC2 = get_rc_data_1D(data_dict, num_sims, rc_diml2, coefl2, num_bins2)

    if avg == 0: # No average to calculate
        pass
    elif avg == 1: # Only average to calculate
        if tot_ene == []:
            raise ValueError('No value list provided to calculate the average value!')
        elif tot_ene != []:
            data_list = tot_ene
            if ene1 != []:
                raise ValueError('Provide two list values to calculate the average '
                                 'value, only one list is needed!')
    elif avg == 2: # About the free energy partitioning
        if tot_ene != [] and ene1 != []:
            if len(tot_ene) == len(ene1):
                #Calculate the exp[beta*(U(ri)-U'(ri))]
                ene2 = [tot_ene[i]-ene1[i] for i in xrange(0, len(tot_ene))]
                ene2min = min(ene2)
                ene2 = [i-ene2min-300.0 for i in ene2]
                data_list = [exp(i/KbT) for i in ene2]
            else:
                raise ValueError('The length of total energy list and paritioning '
                                 'energy list are not the same!')
        elif tot_ene != []:
            raise ValueError('Only total energy list provided, you need to '
                             'provide paritioning energy list as well!')
        elif ene1 != []:
            raise ValueError('Only paritioning energy list provided, you need to '
                             'provide total energy list as well!')

    # Generate the free energy for each bin
    Prob_RC12_N = [[0 for i in xrange(num_binsY)] for j in xrange(num_binsX)]
    Prob_RC12_0 = [[0.0 for i in xrange(num_binsY)] for j in xrange(num_binsX)]
    Prob_RC12_1 = [[0.0 for i in xrange(num_binsY)] for j in xrange(num_binsX)]

    ini_crd = get_string_2D(data_dict, string_seq[0], string_seq[1], rc_diml1,
                            coefl1, rc_diml2, coefl2, figname + '.inistr')
    fin_crd = get_string_2D(data_dict, string_seq[2], string_seq[3], rc_diml1,
                            coefl1, rc_diml2, coefl2, figname + '.laststr')

    count = 0
    kk = 0

    if avg == 0:
        for i in xrange(0, num_sims):
            data_per_sim = len(data_dict[i+1].data)
            for l in xrange(0, data_per_sim):
 
                # For the first dimension
                each_data_RC1 = data_RC1[count]
                each_count_RC1, bins_RCp1 = histogram(each_data_RC1, bins=bins_RC1)
                count_indx1 = list(each_count_RC1).index(1)

                # For the second dimension
                each_data_RC2 = data_RC2[count]
                each_count_RC2, bins_RCp2 = histogram(each_data_RC2, bins=bins_RC2)
                count_indx2 = list(each_count_RC2).index(1)
 
                #Get the biased potentials
                Ubiasl_j = []
                for j in xrange(num_sims):
                    Ubiasl_j.append(Ubiasl[kk])
                    kk = kk + 1
 
                denom = 0.0
                for k in xrange(num_sims):
                    data_per_sim2 = len(data_dict[k+1].data)
                    denom = denom + float(data_per_sim2) * Ubiasl_j[k] * expFx[k]

                Prob_RC12_N[count_indx1][count_indx2] = Prob_RC12_N[count_indx1][count_indx2] + 1
                Prob_RC12_0[count_indx1][count_indx2] = Prob_RC12_0[count_indx1][count_indx2] + 1.0/denom
                count = count + 1

        plot_free_ene_2D(num_binsX, num_binsY, xaxis_RC1, xaxis_RC2, Prob_RC12_N, figname +'.count.pdf',
                         xlabel1, xlabel2, colmap, ini_crd, fin_crd, 0, 1)
        plot_free_ene_2D(num_binsX, num_binsY, xaxis_RC1, xaxis_RC2, Prob_RC12_0, figname,
                         xlabel1, xlabel2, colmap, ini_crd, fin_crd, avg)

    else:
        for i in xrange(0, num_sims):
            data_per_sim = len(data_dict[i+1].data)
            for l in xrange(0, data_per_sim):
 
                # For the first dimension
                each_data_RC1 = data_RC1[count]
                each_count_RC1, bins_RCp1 = histogram(each_data_RC1, bins=bins_RC1)
                count_indx1 = list(each_count_RC1).index(1)
                # For the second dimension
                each_data_RC2 = data_RC2[count]
                each_count_RC2, bins_RCp2 = histogram(each_data_RC2, bins=bins_RC2)
                count_indx2 = list(each_count_RC2).index(1)
 
                #Get the biased potentials
                Ubiasl_j = []
                for j in xrange(num_sims):
                    Ubiasl_j.append(Ubiasl[kk])
                    kk = kk + 1
 
                denom = 0.0
                for k in xrange(num_sims):
                    data_per_sim2 = len(data_dict[k+1].data)
                    denom = denom + float(data_per_sim2) * Ubiasl_j[k] * expFx[k]

                Prob_RC12_N[count_indx1][count_indx2] = Prob_RC12_N[count_indx1][count_indx2] + 1
                Prob_RC12_0[count_indx1][count_indx2] = Prob_RC12_0[count_indx1][count_indx2] + 1.0/denom
                Prob_RC12_1[count_indx1][count_indx2] = Prob_RC12_1[count_indx1][count_indx2] + data_list[count]/denom         
                count = count + 1

        Prob_RC12_2 = [[inf for i in xrange(num_binsY)] for j in xrange(num_binsX)]
        for i in xrange(num_binsY):
            for j in xrange(num_binsX):
                if Prob_RC12_0[j][i] != 0:
                    Prob_RC12_2[j][i] = Prob_RC12_1[j][i]/Prob_RC12_0[j][i]

        plot_free_ene_2D(num_binsX, num_binsY, xaxis_RC1, xaxis_RC2, Prob_RC12_N, figname +'.count.pdf',
                         xlabel1, xlabel2, colmap, ini_crd, fin_crd, 0, 1)
        plot_free_ene_2D(num_binsX, num_binsY, xaxis_RC1, xaxis_RC2, Prob_RC12_2, figname,
                         xlabel1, xlabel2, colmap, ini_crd, fin_crd, avg)

###############################################################################
                       #The N dimension for final string
###############################################################################

def get_string_nD(data_dict, ini_image, fin_image):

    crd = []
    for i in xrange(ini_image, fin_image+1):
        crd.append(data_dict[i].equ_dis)
    crd = array(crd)

    return crd

def get_weights_all(data_dict, num_sims, expFx, Ubiasl):
    """Calculate the weights for each point"""

    weights = []

    kk = 0
    for i in xrange(0, num_sims): #Sum over big N for i
        data_per_sim = len(data_dict[i+1].data)
        for l in xrange(0, data_per_sim): #Sum over little n for l

            # Get the biased potentials
            Ubiasl_j = []
            for j in xrange(num_sims):
                Ubiasl_j.append(Ubiasl[kk])
                kk = kk + 1

            # Obtain the denominator
            denom = 0.0
            for k in xrange(num_sims):
                data_per_sim2 = len(data_dict[k+1].data)
                denom = denom + float(data_per_sim2) * Ubiasl_j[k] * expFx[k]

            weight = 1.0/denom
            weights.append(weight)

    return weights

def gene_free_energy_nD(data_dict, num_sims, string_seq, expFx, Ubiasl, binsize, figname):

    global KbT

    data_all = get_data_all(data_dict, num_sims)
    fin_crd = get_string_nD(data_dict, string_seq[2], string_seq[3])
    weights_all = get_weights_all(data_dict, num_sims, expFx, Ubiasl)

    probs = []
    for i in xrange(len(fin_crd)):
        rangei = [[fin_crd[i][j]-binsize/2.0,fin_crd[i][j]+binsize/2.0] for j in xrange(len(fin_crd[i]))]
        prob_i, edges_i = histogramdd(data_all, range=rangei, weights=weights_all, bins=1)
        prob_i = prob_i.sum()
        probs.append(prob_i)

    n_points = len(fin_crd)
    n_list = range(1, len(fin_crd)+1)
    plot_free_ene_1D(n_points, probs, n_list, figname, 'Image Number')

###############################################################################
                              #The WHAM iteration
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

def wham_iter(num_sims, wham_conv, data_dict, Ubiasl):

    start_time = time.time()

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
        kk=0 #A number to count

        """
        #The code from the WHAM paper
        for k in xrange(0, num_sims):
            ebfk = 0.0
            for i in xrange(0, num_sims):
                data_per_sim = len(data_dict[i+1].data)
                for l in xrange(0, data_per_sim):
                    bottom = 0.0
                    for j in xrange(0, num_sims):
                        bottom = Ubiasl[kk] * expFx_old[j]
                    ebfk = ebfk + Ubiasl[kk]/bottom"""

        for i in xrange(0, num_sims): #Sum over big N for i
            data_per_sim = len(data_dict[i+1].data) 
            for l in xrange(0, data_per_sim): #Sum over little n for l
                denom = 0.0 # Obtain the denominator
                Ubiasl_j = []
                for j in xrange(num_sims):
                   data_per_sim2 = len(data_dict[j+1].data)
                   #denom = denom + Ubiasl[kk] * expFx_old[j]
                   denom = denom + float(data_per_sim2) * Ubiasl[kk] * expFx_old[j]
                   Ubiasl_j.append(Ubiasl[kk])
                   kk = kk + 1
                denom = 1.0/denom
                for k in xrange(num_sims):
                    Fx[k] = Fx[k] + Ubiasl_j[k] * denom

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

    expFx = [exp(i/KbT) for i in Fx]
    write_list('Fx.dat', Fx, num_sims)
    write_list('Fx_prog.dat', Fx_prog, num_sims)

    cost_time = time.time() - start_time
    print("%d Iterations were taken!" %(iter))
    print("It costs %f seconds to finish the WHAM cycle!" %cost_time)

    return expFx

def wham(dirpath, dim, react_paths, first_num_imgs, num_cycles, wham_conv, btstrap=0,
         read_Fx="", read_Ubiasl=""):

    global KbT

    if btstrap not in [0, 1]:
       raise ValueError("bstrap needs to be 0 or 1!")

    # The initial string and final string
    first_i = 1
    first_t = first_num_imgs

    #############################READ THE DATA FILES###########################

    # All of the following list has dimension of iteration * images
    res_crdl = []
    res_forcel = []
    datal = []

    for i in range(1, num_cycles+1):
        window_dir = dirpath + str(i) + '/'
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
            raise ValueError('The number of distances (%d) and (%d) force constants are '
                             'not the same in the step %d!' %(n1, n2, i))

    # Define the final string
    last_num_imgs = n1
    final_t = len(datal)
    final_i = final_t - last_num_imgs + 1
    string_seq = [first_i, first_t, final_i, final_t]

    # Check the three lists are consistent
    if len(res_crdl) == len(res_forcel):
        if len(res_crdl) == len(datal):
            num_sims = len(res_crdl)
        else:
            raise ValueError('The numbers of equlibrium distances and windows are not consistent!')
    else:
        raise ValueError('The number of equlibrium distances and force constants are not consistent!')

    # Apply the fake bootstrapping
    if btstrap == 1:
        datal2 = datal
        seed() # Initialize the random generator
        for i in xrange(num_sims):
            for j in xrange(len(datal[i])):
                k = randint(0, len(datal[i])-1)
                datal2[i][j] = datal[i][k]
        datal = datal2
        del datal2

    # Generate the data_dict, which contains all the data
    data_dict = {}
    for i in range(0, num_sims):
        window_datai = window_data(res_crdl[i], res_forcel[i], datal[i], dim)
        data_dict[i+1] = window_datai

    #Plot the first and last string
    plot_string_1D(first_num_imgs, dim, data_dict, first_i, first_t, react_paths, 'First_string.pdf')
    plot_string_1D(last_num_imgs, dim, data_dict, final_i, final_t, react_paths, 'Last_string.pdf')

    #############################THE WHAM CODE#################################
    #Get the biased potentials

    if read_Ubiasl == "":
        Ubiasl = get_Ubiasl(data_dict, num_sims, dim)
    else:
        Ubiasl = read_list(read_Ubiasl, 1)

    if read_Fx == "":
        expFx = wham_iter(num_sims, wham_conv, data_dict, Ubiasl)
    else:
        Fx = read_dat_file(read_Fx)
	Fx = Fx[0]
        Fx0 = Fx[0] #Normalize the Fx values
        Fx = [Fx[i]-Fx0 for i in xrange(num_sims)]
        expFx = [exp(i/KbT) for i in Fx]

    # RETURN THE DATA FOR FOLLOWING CALCULATIONS
    return data_dict, num_sims, string_seq, expFx, Ubiasl

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
