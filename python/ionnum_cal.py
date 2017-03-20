#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: ionnum_cal.py
from __future__ import print_function
from optparse import OptionParser

parser = OptionParser("Usage: ionnum_cal.py -c protein_charge -d ion_density -x dimension_x -y dimension_y -z dimension_z")
parser.set_defaults(density=0.15)

parser.add_option("-c", dest="charge", type='int',
                  help="Total charge of the protein (without waters or counterions)")
parser.add_option("-d", dest="density", type='float',
                  help="Ion density (Unite: mol/L (or called M))")
parser.add_option("-x", dest="dimx", type='float',
                  help="Dimension along the X axis (Unit: Angstrom)")
parser.add_option("-y", dest="dimy", type='float',
                  help="Dimension along the Y axis (Unit: Angstrom)")
parser.add_option("-z", dest="dimz", type='float',
                  help="Dimension along Z axis (Unit: Angstrom)")
(options, args) = parser.parse_args()

AVOGA_CONSP=6.02214086
volume = options.dimx*options.dimy*options.dimz
counter_num = 0.15 * volume * AVOGA_CONSP * (10.0**-4)
#The 10**-4 Comes from 10**23*10**-27:
#    10**23 from the power of Avogadro's constant
#    10**-27 is the unit transfer between Angstrom**3 and L.
counter_num = round(counter_num, 0)

Na_num = 0
Cl_num = 0
if options.charge > 0:
     Cl_num = counter_num + options.charge
     Na_num = counter_num
elif options.charge < 0:
     Na_num = counter_num - options.charge
     Cl_num = counter_num     
else:
     Na_num = counter_num
     Cl_num = counter_num

print('You need %d Na+ and %d Cl-.' %(Na_num, Cl_num))

