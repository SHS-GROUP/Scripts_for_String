#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: wham_test.py
from __future__ import print_function
from numpy import average, std
from numpy import array, matrix, exp, log, linspace, histogram
import matplotlib.pyplot as plt
import time

start_time = time.time()

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

def write_list(fname, list_name, num_sims):
    writef = open(fname, 'w')
    for i in xrange(1, len(list_name)+1):
        j = list_name[i-1]
        print("%16.7e" %j, end='', file=writef)
        if i%num_sims == 0:
            print(' ', file=writef)
    writef.close()

def read_file(fname):
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

def plot_string(num_imgs, dim, data_dict, ini_wind, final_wind, react_paths, figname):
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
    #plt.show()

def plot_free_ene_1D(x, y, figname, clr_num, rc_label):
    clr = color_dict[clr_num]
    plt.plot(x, y, color=clr, linestyle='-')
    plt.ylabel('Free Energy (kcal/mol)')
    plt.xlabel(rc_label)
    plt.savefig(figname)
    plt.close()

def gene_free_ene_1D(rc_dim, figname, clr_num, rc_label):

    # Get the sampling data for a certain dimension
    data_RC = []

    for i in xrange(0, num_sims):
        data_per_sim = len(data_dict[i+1].data)
        for j in xrange(0, data_per_sim):
            data_RC.append(data_dict[i+1].data[j,rc_dim-1])

    # Determine the plotting bin size
    val_min = min(data_RC) - 0.00001
    val_max = max(data_RC) + 0.00001
    #print(val_min, val_max)

    bins_RC = linspace(val_min, val_max, num_bins+1)
    #print('bins_RC', bins_RC)
    binsize_RC = bins_RC[1] - bins_RC[0]
    plot_bins_RC = bins_RC[0:-1]
    xaxis_RC = [i + binsize_RC/2.0 for i in plot_bins_RC]
    #print('xaxis_RC', xaxis_RC)

    # Generate the free energy for each bin
    Prob_RC = [0.0 for i in xrange(num_bins)]

    for i in xrange(0, num_sims):
        data_per_sim = len(data_dict[i+1].data)
        for j in xrange(0, data_per_sim):
            each_data_RC = data_dict[i+1].data[j,rc_dim-1]
            each_count_RC, bins_RCp = histogram(each_data_RC, bins=bins_RC)
            count_indx = list(each_count_RC).index(1)
            #print(count_indx)

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

    free_ene_RC = [0.0 for ii in xrange(num_bins)]

    for i in xrange(num_bins):
        free_ene_RC[i] = -KbT*log(Prob_RC[i])

    min_free_ene = min(free_ene_RC)
    free_ene_RC = [free_ene - min_free_ene for free_ene  in free_ene_RC]
    plot_free_ene_1D(xaxis_RC, free_ene_RC, figname, clr_num, rc_label)

##############################MAIN PROGRAM#####################################

dim=3 #Dimension of the reaction coordinates
#Reaction path names
react_paths = ['C-H', 'O-H', 'C-O']
num_iter=5 #Number of iteration for the string calculations
num_imgs=18 #Number of images for the string
num_bins=20 #Number of bins
final_i = num_imgs*(num_iter-1) + 1
final_t = final_i + num_imgs - 1

#About the color dictionary
#color_dict = {1: 'b', 2: 'g',
#              3: 'r', 4: 'c',
#              5: 'm', 6: 'y',
#              7: 'y', 8: 'k'}

color_dict = {}
color_disc = [[0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]

if dim <= 7:
    pass
else:
    tmp_num1=(dim - 7)%6
    tmp_num2=(dim - 7)//6
    if tmp_num1 != 0:
        tmp_num2 = tmp_num2 + 1
    for i in xrange(0, len(color_disc1)-1):
       first = color_disc1[i]
       second = color_disc1[i+1]
       for j in xrange(1, tmp_num2):
          color = [first[j] + (second[j] - first[j])* float(j)/float(tmp_num2) for j in xrange(3)]
          color_disc.append(color)
color_disc = sorted(color_disc)

for i in range(0, dim):
    color_dict[i+1] = color_disc[i]

#All of the following list has dimension of iteration * images
res_crdl = []
res_forcel = []
datal = []

for i in range(1, num_iter+1):
    window_dir = './' + str(i) + '/'
    res_crdf = window_dir + 'newconstr.' + str(i) + '.dat'
    res_crd = read_file(res_crdf)
    res_forcef = window_dir + 'force.' + str(i) + '.dat'
    res_force = read_file(res_forcef)
    n1=len(res_crd)
    n2=len(res_force)
    if n1 == n2:
        res_crdl = res_crdl + res_crd
        res_forcel = res_forcel + res_force
        for j in range(1, n1+1):
            datf = window_dir + 'distperframe' + str(j) + '.dat'
            data = read_file(datf)
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
plot_string(num_imgs, dim, data_dict, 1, num_imgs, react_paths, 'First_string.pdf')
plot_string(num_imgs, dim, data_dict, final_i, final_t, react_paths, 'Last_string.pdf')

KbT = 0.596 #Energy unit
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
change = 0.05

while change>0.001:
    # This part is based on the equations 12-15 on the paper
    # of Souaille and Roux on Computer Physics Communications
    # 2001, 135, 40-57
    # Another reference:
    # http://membrane.urmc.rochester.edu/sites/default/files/wham/wham_talk.pdf
    expFx_old = [exp(i/KbT) for i in Fx_old]
    Fx = [0.0 for i in xrange(num_sims)] #Initial free energy
    kk=0
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

cost_time = time.time() - start_time

print("%d Iterations were taken!" %(iter))
print("It costs %f seconds!" %cost_time)

for i in xrange(0, dim):
    plot_dim = i + 1
    clr_num = i + 1
    figname = 'RC' + str(i+1) + '_Free_Energy.pdf'
    gene_free_ene_1D(plot_dim, figname, clr_num, react_paths[i])
