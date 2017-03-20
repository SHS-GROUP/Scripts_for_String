#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
from __future__ import print_function
from optparse import OptionParser
from numpy import zeros, sqrt

parser = OptionParser("Usage: -s step -d dimension -i img -r range --dir running_dir")
parser.add_option("-s", dest="step", type='int',
                  help="Step number")
parser.add_option("-d", dest="dim", type='int',
                  help="Dimension")
parser.add_option("-i", dest="img", type='int',
                  help="Image number")
parser.add_option("-r", dest="rng", type='string',
                  help="Comparison range")
parser.add_option("--dir", dest="rundir", type='string',
                  help="Running directory")
(options, args) = parser.parse_args()

def read_string(fname):
    datm = zeros((options.img, options.dim))
    w_readf = open(fname, 'r')
    count = 1
    for rline in w_readf:
        if count <= options.img:
            line = rline.strip('\n')
            line = line.split()
            for i in xrange(0, options.dim):
                datm[count-1,i] = float(line[i])
        count += 1
    w_readf.close()
    return datm

com_rng = options.rng
com_rng = com_rng.strip(' ')
com_rng = com_rng.split('-')
com_rng = [int(i) for i in com_rng]

# Get the average constraits of the range
avgcons = zeros((options.img, options.dim))
for i in xrange(com_rng[0], com_rng[1]+1):
    dir = options.rundir + '/optstring_calcs' + str(i)
    fname = dir + '/newconstr.dat'
    tempcons = read_string(fname)
    avgcons += tempcons
avgcons = avgcons/float(com_rng[1]-com_rng[0]+1)

# Get the target constraits and calculate the differences
dir = options.rundir + '/optstring_calcs' + str(options.step)
fname = dir + '/newconstr.dat'
tcons = read_string(fname)
diff = avgcons - tcons

# Calculate the RMSD differences along each dimensions and of total
diff = diff**2

sumdiff = sum(diff)
sumsq = sqrt(sumdiff)
print('RMSD along each dimension: ', sumsq)

sumdiff = sum(sum(diff))
sumsq = sqrt(sumdiff)
print('Total RMSD: ', sumsq)

