#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# File name rege_string.py
from __future__ import print_function, division
from numpy import (array, mat, polyfit, polyval, polyder,
                   linspace, sqrt, trapz, zeros)
import matplotlib.pyplot as plt
from string_functions import get_color_dict, get_bpairsl
from optparse import OptionParser

parser = OptionParser("Usage: -i input_file -d dimension --img img_number --nimg new_img_number\n"
                      "       [-f file_suffix] [-o output_file] [-k force_const] [--RST RST_file]\n"
		      "       [--or order]")

parser.set_defaults(f_suffix='dat', outputf= 'newconstr', order=5, rstf='newconstr.RST', fcons=50.0)
parser.add_option("-i", dest="input", type='string',
                  help="Input file name")
parser.add_option("-d", dest="str_dim", type='int',
                  help="String dimension")
parser.add_option("--img", dest="imgs", type='int',
                  help="Number of images")
parser.add_option("--nimg", dest="nimgs", type='int',
                  help="New number of images")
parser.add_option("-f", dest="f_suffix", type='string',
                  help="File name suffix")
parser.add_option("-o", dest="outputf", type='string',
                  help="Output file name prefix")
parser.add_option("-k", dest="fcons", type='float',
                  help="Force constant for AMBER (k in k(r-req)**2) [!ATTENTION, NO 1/2 in the form.]")
parser.add_option("--RST", dest="rstf", type='string',
                  help="RST file name prefix")
#parser.add_option("--prog", dest="program", type='string',
#                  help="program")
parser.add_option("--or", dest="order", type='int',
                  help="Polynomial order for string fitting")
(options, args) = parser.parse_args()

def read_datf(fname, str_dim):
    disl = []
    count = 1
    read_file = open(fname, 'r')
    for rline in read_file:
        if count <= str_dim:
            line = rline.strip('\n')
            line = line.split()
            ati = int(line[0])
            atj = int(line[1])
            dis = float(line[2])
            disl.append(dis)
        count += 1
    read_file.close()
    disl = mat(disl).T
    return disl

bpairsl = get_bpairsl(options.input)

color_dict = get_color_dict(options.str_dim)

avg_val_matrix = [[0.0 for i in range(options.imgs)] for j in range(options.str_dim)]
avg_val_matrix = mat(avg_val_matrix)

for i in xrange(0, options.imgs):
    fn = str(i+1) + '.' + options.f_suffix
    avg_val_matrix[:,i] = read_datf(fn, options.str_dim)

xval = range(1, options.imgs+1)
xval = array(xval)
# X value with high resolution for calculating integral
xval_hr = [0.1*i for i in xrange(10, options.imgs*10+1)] 
xval_hr = array(xval_hr)

# Print for each dimension
pcoefl = []
for dim in xrange(0, options.str_dim):
    yval = avg_val_matrix[dim]
    yval = array(yval)[0]
    pcoef = polyfit(xval, yval, options.order)
    pcoefl.append(pcoef)
    yval_fitted = [polyval(pcoef, i) for i in xval_hr]
    plt.plot(xval, yval, 'ro', xval_hr, yval_fitted, 'b-')
    fig_name = str(dim+1) + '.pdf'
    plt.xticks(xval)
    plt.savefig(fig_name, format='pdf')
    plt.close()

# Print for all the dimensions
for dim in xrange(0, options.str_dim):
    yval = avg_val_matrix[dim]
    yval = array(yval)[0]
    pcoef = polyfit(xval, yval, options.order)
    yval_fitted = [polyval(pcoef, i) for i in xval_hr]
    clr = color_dict[dim+1]
    plt.plot(xval, yval, color=clr, linestyle='', marker='o')
    plt.plot(xval_hr, yval_fitted, 'k-')
plt.xticks(xval)
plt.savefig('total.pdf', format='pdf')
plt.close()

# Calcualte the total length of the string
# Based on the formulation of ds = sqrt((dy1/dx)**2 + (dy2/dx)**2 + ... + (dyn/dx)**2) * dx
# And then string length s = intergral(ds)
#                          = integral (sqrt((dy1/dx)**2 + (dy2/dx)**2 + ... + (dyn/dx)**2) * dx) along x
yval_hr = [0.0 for i in xval_hr]
yval_hr = array(yval_hr)
for dim in xrange(0, options.str_dim):
    yval_hr = yval_hr + (polyval(polyder(pcoefl[dim]), xval_hr))**2
yval_hr_sq = [sqrt(i) for i in yval_hr]
line = trapz(yval_hr_sq, x=xval_hr)

# Redistribute the image points based on the new image numbers
newline = linspace(0, line, options.nimgs)
#print(newline)

# Grid the length along the string
flen = []
for i in xrange(1, len(xval_hr)+1):
    tmp_val1 = trapz(yval_hr_sq[0:i], xval_hr[0:i])
    flen.append(tmp_val1)

# Find corresponding X value index
indl = []
for i in xrange(0, options.nimgs):
    tmp_val2 = newline[i]
    tmp_l1 = [abs(j - tmp_val2) for j in flen]
    minval = min(tmp_l1)
    minind = tmp_l1.index(minval)
    indl.append(minind)

# Assign new values
gmatrx = zeros((options.str_dim, len(xval_hr)))
gmatrx = mat(gmatrx)
newconstr = zeros((options.str_dim, options.nimgs))
newconstr = mat(newconstr)
for i in xrange(0, options.str_dim):
    gmatrx[i,:] = polyval(pcoefl[i], xval_hr)
    newconstr[i,:] = [gmatrx[i,j] for j in indl]

# Print out the new constraints
w_file = open(options.outputf + '.dat', 'w')
for i in xrange(0, options.nimgs):
    for j in xrange(0, options.str_dim):
        print('   %7.4f' %newconstr[j, i],end='', file=w_file)
    print('', file=w_file)
w_file.close()

# Print the new constraints as AMBER input file
w_output = open(options.rstf, 'w')
for i in xrange(0, options.nimgs):
    disl = newconstr[:,i]
    for j in xrange(0, options.str_dim):
        at1 = bpairsl[j][0]
        at2 = bpairsl[j][1]
        print('&rst iat=%d,%d, r1=0., r2=%7.4f, r3=%7.4f, r4=100., rk2=%5.1f, rk3=%5.1f,/'
              %(at1, at2, disl[j], disl[j], options.fcons, options.fcons), file=w_output)
w_output.close()

quit()
