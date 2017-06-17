#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: auto_corr.py
from __future__ import print_function
from numpy import correlate, corrcoef, average, array, allclose, arange
from string_functions import read_dat_file
from matplotlib import pyplot as plt
from optparse import OptionParser

parser = OptionParser("Usage: -f file_name -d dimension")
parser.add_option("-f", dest="fn", type='string',
                  help="Input file name")
parser.add_option("-d", dest="dim", type='int',
                  help="Value dimension")
(options, args) = parser.parse_args()

def estimated_autocorrelation(x):
    """
    http://stackoverflow.com/q/14297012/190597
    http://en.wikipedia.org/wiki/Autocorrelation#Estimation
    """
    n = len(x)
    variance = x.var()
    x = x-x.mean()
    r = correlate(x, x, mode = 'full')[-n:]
    assert allclose(r, array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
    result = r/(variance*(arange(n, 0, -1)))
    return result

def serial_corr(wave, lag=1):
    n = len(wave)
    y1 = wave[lag:]
    y2 = wave[:n-lag]
    corr = corrcoef(y1, y2, ddof=0)[0, 1]
    return corr

def autocorr(wave):
    lags = range(len(wave)//2)
    corrs = [serial_corr(wave, lag) for lag in lags]
    return lags, corrs

#def autocorr(x):
#    result = correlate(x, x, mode='full')
#    return result[result.size/2:]

data = read_dat_file(options.fn)

#print(data)

x = [data[i][options.dim-1] for i in xrange(0, len(data))]

# First Method
xavg = average(x)

x1 = [i-xavg for i in x]

lags, corrs = autocorr(x1)

plt.plot(lags, corrs, 'k-')
plt.savefig(str(options.dim) + '.pdf', format='pdf')
plt.close()

# Second Method
from pandas.tools.plotting import autocorrelation_plot

plt.figure()
autocorrelation_plot(x)
plt.savefig(str(options.dim) + '_pan.pdf', format='pdf')
plt.close()

# Third Method
# http://stackoverflow.com/questions/14297012/estimate-autocorrelation-using-python
x = array(x)
result = estimated_autocorrelation(x)
xval = xrange(1, len(result)+1)
plt.plot(xval, result, 'k-')
plt.savefig(str(options.dim) + '_3.pdf', format='pdf')
plt.close()

