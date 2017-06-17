#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: string_functions.py
from __future__ import print_function
from numpy import average, std, arange
from numpy import array, matrix, exp, log, linspace, histogram
import matplotlib.pyplot as plt
import math
import time

class window_data:
    def __init__(self, equ_dis, constr, data, dim):
        self.equ_dis = equ_dis
        self.constr = constr
        self.data = array(data)
        avg_val_list = []
        std_val_list = []
        for i in xrange(0, dim):
            val = self.data[:,i]
            avg_val = average(val)
            avg_val_list.append(avg_val)
            std_val = std(val)
            std_val_list.append(std_val)
        self.avg = avg_val_list
        self.std = std_val_list

def get_color_dict(dim):
    color_dict = {}
    color_disc = [[0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]

    if dim <= 7:
       pass
    else:
        tmp_num1=(dim - 7)%6
        tmp_num2=(dim - 7)//6 + 1

        if tmp_num1 != 0:
            tmp_num2 = tmp_num2 + 1

        for i in xrange(0, len(color_disc)-1):
           first = color_disc[i]
           second = color_disc[i+1]
           for j in xrange(1, tmp_num2):
              color = [first[k] + (second[k] - first[k]) * float(j)/float(tmp_num2) for k in xrange(3)]
              color_disc.append(color)

    color_disc = sorted(color_disc)

    for i in range(0, dim):
        color_dict[i+1] = color_disc[i]
    return color_dict

def write_list(fname, list_name, num_sims):
    writef = open(fname, 'w')
    for i in xrange(1, len(list_name)+1):
        j = list_name[i-1]
        print("%16.7e" %j, end='', file=writef)
        if i%num_sims == 0:
            print(' ', file=writef)
    writef.close()

def write_xy_lists(fname, xlist, ylist):

    if len(xlist) != len(ylist):
        raise ValueError('The length of xlist and ylist should be the same!')

    writef = open(fname, 'w')
    for i in xrange(len(xlist)):
        j = xlist[i]
        k = ylist[i]
        print("%16.3f %16.3f" %(j, k), file=writef)
    writef.close()

def read_dat_file(fname):
    readf = open(fname, 'r')
    dat_array = []
    ln = 0
    for line in readf:
        line = line.strip('\n')
        line = line.split()
        line = [float(i) for i in line]
        if len(line) != 0:
            dat_array.append(line)
    readf.close()
    return dat_array

def get_rc_data_1D(data_dict, num_sims, rc_diml, coefl, num_bins):

    if len(rc_diml) != len(coefl):
        raise ValueError("The length of the coordinate dimensions and "
                         "coefficients in the linear combination of "
                         "reaction coordinate are not equal!")

    # Get the sampling data for a certain dimension
    data_RC = []
    data_RCl = []
    for i in xrange(0, num_sims):
        data_per_sim = len(data_dict[i+1].data)
        data_RC_per_sim = []
        for j in xrange(0, data_per_sim):
            data_val = 0.0
            for k in xrange(len(rc_diml)):
                tmp_data = data_dict[i+1].data[j,rc_diml[k]-1]
                tmp_data = tmp_data * float(coefl[k])
                data_val = data_val + tmp_data
            data_RC.append(data_val)
            data_RC_per_sim.append(data_val)
        data_RCl.append(data_RC_per_sim)

    # Determine the plotting bin size
    val_min = min(data_RC) - 0.00001
    val_max = max(data_RC) + 0.00001

    bins_RC = linspace(val_min, val_max, num_bins+1)
    binsize_RC = bins_RC[1] - bins_RC[0]
    plot_bins_RC = bins_RC[0:-1]
    xaxis_RC = [i + binsize_RC/2.0 for i in plot_bins_RC]

    # Global binning
    color_dict = get_color_dict(num_sims)

    for i in xrange(0, num_sims):
        data_RC_per_sim = data_RCl[i]
        count_RC, bins_RCp = histogram(data_RC_per_sim, bins=bins_RC)
        clr = color_dict[i+1]
        plt.plot(xaxis_RC, count_RC, color=clr, linestyle='-')

    figname = ''
    for i in xrange(0, len(rc_diml)):
        figname = figname + str(round(coefl[i],1)) + 'R' + str(rc_diml[i]) 
    figname = figname + 'globin.pdf'
    plt.savefig(figname)
    plt.close()

    return data_RC, bins_RC, xaxis_RC

def cal_dis(crds, ati, atj):
    sq = 0.0
    for i in xrange(3):
        crd1 = crds[ati-1][i]
        crd2 = crds[atj-1][i]
        sq += (crd1 - crd2)**2
    dis = math.sqrt(sq)
    return dis

def get_bpairsl(fname):
    bpairsl = []
    inputf = open(fname, 'r')
    for line in inputf:
        line = line.split()
        if '\n' in line:
            line.remove('\n')
        if ',' in line:
            line.remove(',')
        if '' in line:
            line.remove('')
        if ' ' in line:
            line.remove(' ')
        if ':' in line:
            line.remove(':')
        # Blank line
        if (len(line) == 0):
            continue
        # Comment
        elif (line[0][0] == '#'):
            continue
        elif len(line) == 1:
            bpairs0 = line
            for i in bpairs0:
                bpairs = i.split('-')
                try:
                    bpairs = [int(i) for i in bpairs]
                except:
                    raise ValueError('Should be integer numbers in the '
                                   'two ends of dash symbol.')
                if len(bpairs) != 2:
                    raise ValueError('Should be only two numbers in the pairs!')
                bpairs = tuple(bpairs)
                bpairsl.append(bpairs)
    inputf.close()
    return bpairsl

def read_list(fname, row):
    data_list = []
    r_file = open(fname, 'r')
    for rline in r_file:
        line = rline.strip('\n')
        line = line.strip(' ')
        line = line.split()
        data_list.append(float(line[row-1]))
    r_file.close()
    return data_list

