#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: cal_dis.py
from __future__ import print_function
import MDAnalysis
from numpy.linalg import norm
from numpy import array, dot
import parmed as pmd
from string_functions import get_bpairsl, cal_dis
from optparse import OptionParser

parser = OptionParser("Usage:  -i input_file -p top_file -c traj_file -o output_file")

parser.add_option("-i", dest="input", type='string',
                  help="Input file name")
parser.add_option("-p", dest="topf", type='string',
                  help="Topology file name")
parser.add_option("-c", dest="trajf", type='string',
                  help="Trajectory file name")
parser.add_option("-o", dest="output", type='string',
                  help="Output file name")
(options, args) = parser.parse_args()

###############################################################################
# Read the input file
###############################################################################

exlist = []
scalist = []
point = 1
donor = 1
acceptor = 1
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
    #---------------------------Required Variables-----------------------------
    #exclude list
    elif line[0].lower() == 'exclude':
        if len(line) >= 2:
            try:
                exlist = line[1:]
                exlist = [int(i) for i in exlist]
            except:
                raise ValueError('exclude need to be integer number(s).')
    #scale list
    elif line[0].lower() == 'scale':
        if len(line) >= 2:
            try:
                scalist = line[1:]
                scalist = [int(i) for i in scalist]
            except:
                raise ValueError('scale need to be integer number(s).')
    #point list
    elif line[0].lower() == 'point':
        if len(line) == 2:
            try:
                point = line[1]
                point = int(point)
            except:
                raise ValueError('point need to be integer number(s).')
    #donor
    elif line[0].lower() == 'donor':
        if len(line) == 2:
            try:
                donor = line[1]
                donor = int(donor)
            except:
                raise ValueError('donor need to be integer number(s).')
    #acceptor
    elif line[0].lower() == 'acceptor':
        if len(line) == 2:
            try:
                acceptor = line[1]
                acceptor = int(acceptor)
            except:
                raise ValueError('acceptor need to be integer number(s).')
    #bpairsl
    elif line[0].lower() == 'bpairs':
        if len(line) >= 2:
            addbpairs0 = line[1:]
            for i in addbpairs0:
                pairs = i.split('-')
                try:
                    pairs = [int(i) for i in pairs]
                except:
                    raise ValueError('Should be integer numbers in the '
                                      'two ends of dash symbol.')
                if len(pairs) != 2:
                    raise ValueError('Should be only two numbers in the pairs!')
                pairs = tuple(pairs)
                bpairsl.append(pairs)
        else:
            raise ValueError('bpairs need to be provided.')

print("exclude=", exlist)
print("scale=", scalist)
print("point=", point)
print("donor=", donor)
print("acceptor=", acceptor)
print("bpairs=", bpairsl)

###############################################################################
# Read the CHELPG charges from the output file
###############################################################################

num_atms = 122
hd1 = range(7786, 7797)
hd1.append(7799)
hd2 = range(7875, 7886)
hd2.append(7888)
hd3 = range(10889, 10900)
hd3.append(10902)
an1 = range(10943, 10951)
an1.append(10953)
rest = range(13257, 13334)
tot_list = hd1 + hd2 + hd3 + an1 + rest
mschg_dict = {}
for i in xrange(1, len(tot_list)+1):
    mschg_dict[i] = tot_list[i-1]

def read_chg(fname, num_atms, mol_top, mschg_dict):

    ln = 1
    hasesp1 = 0
    fp = open(fname, 'r')
    for line in fp:
        if 'Ground-State ChElPG Net Atomic Charges' in line:
            hasesp1 = hasesp1 + 1
            bln = ln + 4
        ln = ln + 1
    fp.close()

    if hasesp1 > 0:
        pass
    else:
       raise pymsmtError('There is no \'Ground-State ChElPG Net Atomic Charges\''
             ' found in the Q-Chem output file. '
             'Please check whether the Q-Chem jobs are '
             'finished normally, and whether you are using the '
             'correct output file.')

    chg_dict = {}
    ln = 1
    fp = open(fname, 'r')
    for line0 in fp:
        if (ln >= bln) and (ln <= bln+num_atms-1):
            line = line0.split()
            chg_dict[int(line[0])] = float(line[2])
        ln = ln + 1
    fp.close()

    #print(len(mschg_dict))
    #print(len(chg_dict))
    #print(num_atms)

    if len(mschg_dict) != len(chg_dict) or len(mschg_dict) != num_atms:
        raise ValueError('The length of charge dictionaries are different!')

    for i in xrange(1, num_atms+1):
        mol_top.atoms[mschg_dict[i]-1]._charge = chg_dict[i]

    return mol_top

###############################################################################
# Generate the scaling dictionary
###############################################################################

ke = 8.9875518 #* 10^9 N*m^2*C^-2
e = 1.6021765 #* 10^-19 C
constant = ke * e # 10^9 * 10^-19 * 10^20 = 10^10 N/C = 10^10 V/m = 10^4 MV/m = 10^2 MV/cm

def cal_efield(crdx, crdy, chgx, scale):
    #print(crdx, crdy, chgx, scale)
    global constant
    crdr = crdy - crdx
    #crdr = crdx - crdy
    r = norm(crdr)
    r2 = r * r
    Felec = chgx * constant / r2
    arrow = Felec * scale * crdr/r #use crdr/r to transfer it to a vector
    #print(arrow)
    return arrow

# Read the topology and trajectory files
mol_top = pmd.load_file(options.topf)
mol = MDAnalysis.Universe(options.topf, options.trajf)

# Generate the scaling dict for the electrostatic
scale_dict = {}
for i in xrange(1, len(mol_top.atoms)+1):
    if i in exlist:
        scale_dict[i] = 0.00
    elif i in scalist:
        scale_dict[i] = 1.00/1.20
    else:
        scale_dict[i] = 1.00

###############################################################################
# Calculate the electrostatic forces
###############################################################################

hd1 = range(7782, 7800) #Total HID atoms
hd2 = range(7871, 7889) #Total HID atoms
hd3 = range(10885, 10903) #Total HID atoms
an1 = range(10939, 10954) #Total ASN atoms
ie1 = range(13257, 13280) #Total ILE atoms + COH in SER838
feo = range(13280, 13283) #Fe+OH
lid = range(13283, 13334) #Total LID
tot_num = hd1 + hd2 + hd3 + an1 + ie1 + feo + lid

#ms_atids = [7782, 7783, 7784, 7785, 7786, 7787, 7788, 7789, 7790, 7791, 7792, 7793, 7794, 7795, 7796, 7797, 7798, 7799, 7871, 7872, 7873, 7874, 7875, 7876, 7877, 7878, 7879, 7880, 7881, 7882, 7883, 7884, 7885, 7886, 7887, 7888, 10885, 10886, 10887, 10888, 10889, 10890, 10891, 10892, 10893, 10894, 10895, 10896, 10897, 10898, 10899, 10900, 10901, 10902, 10939, 10940, 10941, 10942, 10943, 10944, 10945, 10946, 10947, 10948, 10949, 10950, 10951, 10952, 10953, 13257, 13258, 13259, 13260, 13261, 13262, 13263, 13264, 13265, 13266, 13267, 13268, 13269, 13270, 13271, 13272, 13273, 13274, 13275, 13276, 13277, 13278, 13279, 13280, 13281, 13282, 13283, 13284, 13285, 13286, 13287, 13288, 13289, 13290, 13291, 13292, 13293, 13294, 13295, 13296, 13297, 13298, 13299, 13300, 13301, 13302, 13303, 13304, 13305, 13306, 13307, 13308, 13309, 13310, 13311, 13312, 13313, 13314, 13315, 13316, 13317, 13318, 13319, 13320, 13321, 13322, 13323, 13324, 13325, 13326, 13327, 13328, 13329, 13330, 13331, 13332, 13333]

# Calculate the distance and electric field strength
w_output = open(options.output, 'w')
count = 1

for traj in mol.trajectory:
    arrow0 = array([0.0, 0.0, 0.0])
    arrowqm = array([0.0, 0.0, 0.0])
    arrowmm = array([0.0, 0.0, 0.0])
    arrow_hd1 = array([0.0, 0.0, 0.0])
    arrow_hd2 = array([0.0, 0.0, 0.0])
    arrow_hd3 = array([0.0, 0.0, 0.0])
    arrow_an1 = array([0.0, 0.0, 0.0])
    arrow_ie1 = array([0.0, 0.0, 0.0])
    arrow_feo = array([0.0, 0.0, 0.0])
    arrow_lid = array([0.0, 0.0, 0.0])

    fname = './saved_outputs/qchem.out_' + str(count)
    mol_top2 = read_chg(fname, num_atms, mol_top, mschg_dict)

    dcrd = mol.atoms[donor-1].position
    acrd = mol.atoms[acceptor-1].position
    damd = (dcrd + acrd) * 0.5
    davec = acrd - dcrd
    davec2 = davec / norm(davec)

    for i in xrange(1, len(mol_top.atoms)+1):
        if i != point:
            crdx = mol.atoms[i-1].position
            crdy = damd
            chgx = mol_top2.atoms[i-1]._charge
            scale = scale_dict[i]
            arrow1 = cal_efield(crdx, crdy, chgx, scale)

            arrow0 += arrow1

            if i in tot_num:
                arrowqm += arrow1
                if i in hd1:
                    arrow_hd1 += arrow1
                elif i in hd2:
                    arrow_hd2 += arrow1
                elif i in hd3:
                    arrow_hd3 += arrow1
                elif i in an1:
                    arrow_an1 += arrow1
                elif i in ie1:
                    arrow_ie1 += arrow1
                elif i in feo:
                    arrow_feo += arrow1
                elif i in lid:
                    arrow_lid += arrow1
            else:
                arrowmm += arrow1

    for bpairs in bpairsl:
        ati = bpairs[0]
        atj = bpairs[1]
        at1 = mol.atoms[ati-1]
        at2 = mol.atoms[atj-1]
        r = at1.position - at2.position
        dis = norm(r)
        print(' %7.4f' %dis, end='', file=w_output)

    arrow0 = 100.0 * dot(arrow0, davec2)
    arrowqm = 100.0 * dot(arrowqm, davec2)
    arrowmm = 100.0 * dot(arrowmm, davec2)
    arrow_hd1 = 100.0 * dot(arrow_hd1, davec2)
    arrow_hd2 = 100.0 * dot(arrow_hd2, davec2)
    arrow_hd3 = 100.0 * dot(arrow_hd3, davec2)
    arrow_an1 = 100.0 * dot(arrow_an1, davec2)
    arrow_ie1 = 100.0 * dot(arrow_ie1, davec2)
    arrow_feo = 100.0 * dot(arrow_feo, davec2)
    arrow_lid = 100.0 * dot(arrow_lid, davec2)

    print(' %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e'  %(arrow0, arrowqm, arrowmm, arrow_hd1, arrow_hd2, arrow_hd3, arrow_an1, arrow_ie1, arrow_feo, arrow_lid), file=w_output)
w_output.close()

quit()

