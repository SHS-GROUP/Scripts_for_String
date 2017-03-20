#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: cal_dis.py
from __future__ import print_function
import MDAnalysis
from numpy.linalg import norm
import parmed as pmd
import math
from optparse import OptionParser

parser = OptionParser("Usage:  -i input_file -p top_file -c traj_file -o output_file -t total_distance_file [-k force_const] [--RST RST_file_name] \n")
#                      "        [-d decimal_places]")
#parser.set_defaults(decimal=3)
parser.set_defaults(rstf='NONONO', fcons=50.0)
parser.add_option("-i", dest="input", type='string',
                  help="Input file name")
parser.add_option("-p", dest="topf", type='string',
                  help="Topology file name")
parser.add_option("-c", dest="trajf", type='string',
                  help="Trajectory file name")
parser.add_option("-o", dest="output", type='string',
                  help="Output file name")
parser.add_option("--RST", dest="rstf", type='string',
                  help="RST file name")
parser.add_option("-k", dest="fcons", type='float',
                  help="Force constant for AMBER (k in k(r-req)**2) [!ATTENTION, NO 1/2 in the form.]")
parser.add_option("-t", dest="output2", type='string',
                  help="Output file name 2, output all the distances.")
#parser.add_option("-d", dest="decimal", type='int',
#                  help="Decimal places keep")
#parser.add_option("--prog", dest="program", type='string',
#                  help="Program: amber/charmm")
(options, args) = parser.parse_args()

def cal_dis(crds, ati, atj):
    sq = 0.0
    for i in xrange(3):
        crd1 = crds[ati-1][i]
        crd2 = crds[atj-1][i]
        sq += (crd1 - crd2)**2
    dis = math.sqrt(sq)
    return dis

bpairsl = []
inputf = open(options.input, 'r')
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

if options.trajf[-4:] == '.dcd': # If it is the CHARMM .dcd binary file
    w_output = open(options.output, 'w')

    if options.rstf != 'NONONO':
        w_rstf = open(options.rstf, 'w')

    mol = MDAnalysis.Universe(options.topf, options.trajf)
    print('There is/are %d frame(s) considered.' %mol.trajectory.n_frames)
    for bpairs in bpairsl:
        ati = bpairs[0]
        atj = bpairs[1]
        disl = []
        for traj in mol.trajectory:
            at1 = mol.atoms[ati-1]
            at2 = mol.atoms[atj-1]
            r = at1.position - at2.position
            dis = norm(r)
            disl.append(dis)
        avgdis = sum(disl)/float(len(disl))
        print('(Average) Distance between Atom %s(ID %d) and Atom %s(ID %d): %7.4f' %(at1.name, ati, at2.name, atj, avgdis))
        print('%d %d %7.4f' %(ati, atj, avgdis), file=w_output)
        #if options.program == 'amber':
        #    print('&rst iat=%d,%d, r1=0., r2=%6.3f, r3=%6.3f, r4=100., rk2=100.0, rk3=100.0,/'
        #          %(ati, atj, avgdis, avgdis), file=w_output)
        #elif options.program == 'charmm':
        #    print('RESDistance KVAL 100 RVAL %6.3f 1.0 bynum %d bynum %d'
        #          %(avgdis, ati, atj), file=w_output)

        if options.rstf != 'NONONO':
                print('&rst iat=%d,%d, r1=0., r2=%7.4f, r3=%7.4f, r4=100., rk2=%5.1f, rk3=%5.1f,/'
                %(ati, atj, avgdis, avgdis, options.fcons, options.fcons), file=w_rstf)

    if options.rstf != 'NONONO':
        w_rstf.close()

    w_output.close()

    if options.output2 is not None:
        w_output2 = open(options.output2, 'w')
        for traj in mol.trajectory:
            for bpairs in bpairsl:
                ati = bpairs[0]
                atj = bpairs[1]
                at1 = mol.atoms[ati-1]
                at2 = mol.atoms[atj-1]
                r = at1.position - at2.position
                dis = norm(r)
                print(' %7.4f' %dis, end='', file=w_output2)
            print('', file=w_output2)
        w_output2.close()

else: # Other formats such as .inpcrd .rst .netcdf .crd .restart
    if options.topf[-4:] == '.prmtop':
    	mol_top = pmd.load_file(options.topf)
    elif options.topf[-4:] == '.psf':
	mol_top = pmd.load_file(options.trajf)

    mol_crd = pmd.load_file(options.trajf)
    crdsl = mol_crd.coordinates
    print('There is/are %d frame(s) considered.' %len(crdsl))
    w_output = open(options.output, 'w')

    if options.rstf != 'NONONO':
        w_rstf = open(options.rstf, 'w')

    for bpairs in bpairsl:
        ati = bpairs[0]
        atj = bpairs[1]
        disl = []
        for crds in crdsl:
            dis = cal_dis(crds, ati, atj)
            disl.append(dis)
        avgdis = sum(disl)/float(len(disl))
	if options.topf[-4:] == '.prmtop':
            atiname = mol_top.atoms[ati-1].name
            atjname = mol_top.atoms[atj-1].name
        elif options.topf[-4:] == '.psf':
            atiname = mol_top.atname[ati-1]
            atjname = mol_top.atname[atj-1]

        print('(Average) Distance between Atom %s(ID %d) and Atom %s(ID %d): %7.4f' %(atiname, ati, atjname, atj, avgdis))
        print('%d %d %7.4f' %(ati, atj, avgdis), file=w_output)
        #if options.program == 'amber':
        #    print('&rst iat=%d,%d, r1=0., r2=%6.3f, r3=%6.3f, r4=100., rk2=100.0, rk3=100.0,/'
        #          %(ati, atj, avgdis, avgdis), file=w_output)
        #elif options.program == 'charmm':
        #    print('RESDistance KVAL 100 RVAL %6.3f 1.0 bynum %d bynum %d'
        #          %(avgdis, ati, atj), file=w_output)

	if options.rstf != 'NONONO':
        	print('&rst iat=%d,%d, r1=0., r2=%7.4f, r3=%7.4f, r4=100., rk2=%5.1f, rk3=%5.1f,/'
               	%(ati, atj, avgdis, avgdis, options.fcons, options.fcons), file=w_rstf)

    if options.rstf != 'NONONO':
	w_rstf.close()

    w_output.close()

    if options.output2 is not None:
        w_output2 = open(options.output2, 'w')
        for crds in crdsl:       
            for bpairs in bpairsl:
                ati = bpairs[0]
                atj = bpairs[1]
                dis = cal_dis(crds, ati, atj)
                print(' %7.4f' %dis, end='', file=w_output2)
            print('', file=w_output2)
        w_output2.close()

