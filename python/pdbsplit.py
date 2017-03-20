#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: pdbsplit.py
from __future__ import print_function
from msmtmol.readpdb import get_atominfo_fpdb
from optparse import OptionParser
import os

parser = OptionParser("Usage: pdbsplit.py -i input_file -p pdb_file --rtf RTF_file --prm PRM_file \n"
                      "                  [-w water] [--renum renum]")

parser.set_defaults(water=0,renum=0)

parser.add_option("-i", dest="inputf", type='string',
                  help="Input file name")
parser.add_option("-p", dest="pdbf", type='string',
                  help="Input PDB file name")
parser.add_option("--rtf", dest="rtf", type='string',
                  help="The RTF file address")
parser.add_option("--prm", dest="prm", type='string',
                  help="The PRM file address")
parser.add_option("-w", dest="water", type='int',
                  help="Maximum water in each water file (0 means not "
                       "automatically split PDB files)")
parser.add_option("--renum", dest="renum", type='int',
                  help="Use pdb4amber to renumber each PDB file or not (0/1)")
#parser.add_option("--charmm", dest="charmm", type='int',
#                  help="Generate CHARMM input files or not (0/1)")

(options, args) = parser.parse_args()

def write_charmm_inp(pdbf, charmmf, rtf, prm):
    resname = pdbf.rstrip('.pdb')
    w_charmmf = open(charmmf, 'w')
    print("* add H to PDB file", file=w_charmmf)
    print("*", file=w_charmmf)
    print("prnlev 5", file=w_charmmf)
    print("wrnlev 8", file=w_charmmf)
    print("!bomlev -2", file=w_charmmf)
    print("", file=w_charmmf)
    print("! read topology and parameter files", file=w_charmmf)
    print("", file=w_charmmf)
    print("read rtf card name \"%s\" unit 1" %rtf, file=w_charmmf)
    print("read para card name \"%s\" unit 1" %prm, file=w_charmmf)
    print("", file=w_charmmf)
    print("open unit 1 card read name \"%s\"" %pdbf, file=w_charmmf)
    print("read sequ pdb unit 1", file=w_charmmf)
    print("close unit 1", file=w_charmmf)
    print("gener %s first none last none" %resname, file=w_charmmf)
    print("", file=w_charmmf)
    print("! Read protein coord from the PDB coordinate file", file=w_charmmf)
    print("open unit 1 card read name \"%s\"" %pdbf, file=w_charmmf)
    print("read coor pdb unit 1", file=w_charmmf)
    print("close unit 1", file=w_charmmf)
    print("", file=w_charmmf)
    print("print coor", file=w_charmmf)
    print("", file=w_charmmf)
    print("ic para", file=w_charmmf)
    print("ic fill preserve", file=w_charmmf)
    print("ic build", file=w_charmmf)
    print("!hbuild sele all end", file=w_charmmf)
    print("", file=w_charmmf)
    print("!coor print", file=w_charmmf)
    print("", file=w_charmmf)
    print("! write out the protein structure file (psf) and", file=w_charmmf)
    print("! the coordinate file in pdb and crd format.", file=w_charmmf)
    print("", file=w_charmmf)
    print("write psf card name %s.psf unit 1" %resname, file=w_charmmf)
    print("", file=w_charmmf)
    print("write coor pdb name %s_new.pdb unit 1" %resname, file=w_charmmf)
    print("", file=w_charmmf)
    print("write coor card name %s.crd unit 1" %resname, file=w_charmmf)
    print("", file=w_charmmf)
    print("stop", file=w_charmmf)
    print("", file=w_charmmf)
    print("END", file=w_charmmf)
    w_charmmf.close()

def write_charmm_combf(bpairsl, charmmf, rtf, prm):
    finalf = charmmf.rstrip('.in')
    w_charmmf = open(charmmf, 'w')
    print("* add H to PDB file", file=w_charmmf)
    print("*", file=w_charmmf)
    print("prnlev 5", file=w_charmmf)
    print("wrnlev 8", file=w_charmmf)
    print("!bomlev -2", file=w_charmmf)
    print("", file=w_charmmf)
    print("! read topology and parameter files", file=w_charmmf)
    print("", file=w_charmmf)
    print("read rtf card name \"%s\" unit 1" %rtf, file=w_charmmf)
    print("read para card name \"%s\" unit 1" %prm, file=w_charmmf)
    print("", file=w_charmmf)
    fname = bpairsl[0][2]
    resname = fname.rstrip('.pdb')
    print("read psf card name %s.psf unit 1" %resname, file=w_charmmf)
    print("read coor card name %s.crd unit 1" %resname, file=w_charmmf)
    print("", file=w_charmmf)
    for i in xrange(1, len(bpairsl)):
        fname = bpairsl[i][2]
        resname = fname.rstrip('.pdb')
        print("read psf card name %s.psf unit 1 append" %resname, file=w_charmmf)
        print("read coor card name %s.crd unit 1 append" %resname, file=w_charmmf)
        print("", file=w_charmmf)
    print("write psf card name %s.psf unit 1" %finalf, file=w_charmmf)
    print("write coor card name %s.crd unit 1" %finalf, file=w_charmmf)
    print("", file=w_charmmf)
    print("stop", file=w_charmmf)
    print("", file=w_charmmf)
    print("END", file=w_charmmf)

ATOM_DICT = {'H5\'1' : 'H5\'',
             'H5\'2' : 'H5\'\'',
             'H2\'1' : 'H2\'\'',
             'HO\'2' : 'H2\'',
             '\'HO3' : 'H3T',
             'HO5\'' : 'H5T',
             'OXT'   : 'OT2',
            }

RES_DICT = {
	'U'  : 'URA',
	'A'  : 'ADE',
	'C'  : 'CYT',
	'G'  : 'GUA',
	'DA' : 'ADE',
	'DG' : 'GUA',
	'DC' : 'CYT',
	'DT' : 'THY',
	'DT3': 'THY',
	'DT5': 'THY'
           }

RES_ATOM_DICT = {
	('ALA', 'H') : ('ALA', 'HN'),
	('ARG', 'H') : ('ARG', 'HN'),
	('ARG', 'HB2') : ('ARG', 'HB1'),
        ('ARG', 'HB3') : ('ARG', 'HB2'),
	('ARG', 'HG2') : ('ARG', 'HG1'),
        ('ARG', 'HG3') : ('ARG', 'HG2'),
        ('ARG', 'HD2') : ('ARG', 'HD1'),
        ('ARG', 'HD3') : ('ARG', 'HD2'),
     	('ASN', 'H')  : ('ASN', 'HN'),
	('ASN', 'HB2') : ('ASN', 'HB1'),
        ('ASN', 'HB3') : ('ASN', 'HB2'),
	('ASP', 'H') : ('ASP', 'HN'),
	('ASP', 'HB2') : ('ASP', 'HB1'),
        ('ASP', 'HB3') : ('ASP', 'HB2'),
	('CYM', 'H') : ('CYM', 'HN'),
	('CYM', 'HB2') : ('CYM', 'HB1'),
	('CYM', 'HB3') : ('CYM', 'HB2'),
	('CYS', 'H') : ('CYS', 'HN'),
        ('CYS', 'HB2') : ('CYS', 'HB1'),
        ('CYS', 'HB3') : ('CYS', 'HB2'),
	('CYS', 'HG') : ('CYS', 'HG1'),
	('CYX', 'H') : ('CYX', 'HN'),
        ('CYX', 'HB2') : ('CYX', 'HB1'),
        ('CYX', 'HB3') : ('CYX', 'HB2'),
	('GLN', 'H') : ('GLN', 'HN'),
	('GLN', 'HB2') : ('GLN', 'HB1'),
	('GLN', 'HB3') : ('GLN', 'HB2'),
	('GLN', 'HG2') : ('GLN', 'HG1'),
	('GLN', 'HG3') : ('GLN', 'HG2'),
	('GLU', 'H') : ('GLU', 'HN'),
        ('GLU', 'HB2') : ('GLU', 'HB1'),
        ('GLU', 'HB3') : ('GLU', 'HB2'),
        ('GLU', 'HG2') : ('GLU', 'HG1'),
        ('GLU', 'HG3') : ('GLU', 'HG2'),
	('GLY', 'H') : ('GLY', 'HN'),
	('GLY', 'HA2') : ('GLY', 'HA1'),
        ('GLY', 'HA3') : ('GLY', 'HA2'),
	('HID', 'H') : ('HID', 'HN'),
	('HID', 'HB2') : ('HID', 'HB1'),
        ('HID', 'HB3') : ('HID', 'HB2'),
	('HIE', 'H') : ('HIE', 'HN'),
	('HIE', 'HB2') : ('HIE', 'HB1'),
        ('HIE', 'HB3') : ('HIE', 'HB2'),
	('HIP', 'H') : ('HIP', 'HN'),
	('HIP', 'HB2') : ('HIP', 'HB1'),
        ('HIP', 'HB3') : ('HIP', 'HB2'),
	('ILE', 'H') : ('ILE', 'HN'),
	('ILE', 'CD1') : ('ILE', 'CD'),
	('ILE', 'HD11') : ('ILE', 'HD1'),
	('ILE', 'HD12') : ('ILE', 'HD2'),
	('ILE', 'HD13') : ('ILE', 'HD3'),
	('ILE', 'HG12') : ('ILE', 'HG11'),
        ('ILE', 'HG13') : ('ILE', 'HG12'),
        ('LEU', 'H') : ('LEU', 'HN'),
        ('LEU', 'HB2') : ('LEU', 'HB1'),
        ('LEU', 'HB3') : ('LEU', 'HB2'),
	('LYS', 'H') : ('LYS', 'HN'),
	('LYS', 'HB2') : ('LYS', 'HB1'),
        ('LYS', 'HB3') : ('LYS', 'HB2'),
	('LYS', 'HG2') : ('LYS', 'HG1'),
        ('LYS', 'HG3') : ('LYS', 'HG2'),
	('LYS', 'HD2') : ('LYS', 'HD1'),
        ('LYS', 'HD3') : ('LYS', 'HD2'),
	('LYS', 'HE2') : ('LYS', 'HE1'),
        ('LYS', 'HE3') : ('LYS', 'HE2'),
	('MET', 'H') : ('MET', 'HN'),
       	('MET', 'HB2') : ('MET', 'HB1'),
        ('MET', 'HB3') : ('MET', 'HB2'),
	('MET', 'HG2') : ('MET', 'HG1'),
        ('MET', 'HG3') : ('MET', 'HG2'),
	('PHE', 'H') : ('PHE', 'HN'),
        ('PHE', 'HB2') : ('PHE', 'HB1'),
        ('PHE', 'HB3') : ('PHE', 'HB2'),
	('PRO', 'HD2') : ('PRO', 'HD1'),
        ('PRO', 'HD3') : ('PRO', 'HD2'),
	('PRO', 'HG2') : ('PRO', 'HG1'),
        ('PRO', 'HG3') : ('PRO', 'HG2'),
	('PRO', 'HB2') : ('PRO', 'HB1'),
        ('PRO', 'HB3') : ('PRO', 'HB2'),
	('SER', 'H') : ('SER', 'HN'),
	('SER', 'HB2') : ('SER', 'HB1'),
        ('SER', 'HB3') : ('SER', 'HB2'),
	('SER', 'HG') : ('SER', 'HG1'),
	('THR', 'H') : ('THR', 'HN'),
	('TRP', 'H') : ('THR', 'HN'),
	('TRP', 'HB2') : ('TRP', 'HB1'),
        ('TRP', 'HB3') : ('TRP', 'HB2'),
	('TYR', 'H') : ('TYR', 'HN'),
	('TYR', 'HB2') : ('TYR', 'HB1'),
        ('TYR', 'HB3') : ('TYR', 'HB2'),
	('VAL', 'H') : ('VAL', 'HN'),
        ('WAT', 'O')   : ('TIP3', 'OH2'),
        ('WAT', 'H1')  : ('TIP3', 'H1'),
        ('WAT', 'H2')  : ('TIP3', 'H2'),
	('Na+', 'Na+') : ('SOD', 'SOD'),
	('Cl-', 'Cl-') : ('CLA', 'CLA'),
		}

"""
RES_ATOM_DICT = {
        ('WAT', 'O')   : ('TIP3', 'OH2'),
        ('WAT', 'H1')  : ('TIP3', 'H1'),
        ('WAT', 'H2')  : ('TIP3', 'H2')
                }
"""

def charmm_mol(mol, atids):
    global ATOM_DICT
    global RES_ATOM_DICT
    for i in atids:
        atname = mol.atoms[i].atname
        resname = mol.atoms[i].resname
        if atname in ATOM_DICT.keys():
            mol.atoms[i].atname = ATOM_DICT[atname]
        elif (resname, atname) in RES_ATOM_DICT.keys():
            mol.atoms[i].resname = RES_ATOM_DICT[(resname, atname)][0]
            mol.atoms[i].atname = RES_ATOM_DICT[(resname, atname)][1]
    return mol

bpairsl = []
inputf = open(options.inputf, 'r')
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
    elif len(line) == 2:
        bpairs0 = line[0]
        bpairs = bpairs0.split('-')
        try:
            bpairs = [int(i) for i in bpairs]
        except:
            raise ValueError('Should be integer numbers in the '
                           'two ends of dash symbol.')
        if len(bpairs) != 2:
            raise ValueError('Should be only two numbers in the pairs!')
        if bpairs[0] > bpairs[1]:
            raise ValueError('The later number should be bigger than the '
                'former number in the range specifications!')

        bpairs = [bpairs[0], bpairs[1], line[1]]
        bpairsl.append(bpairs)
print(bpairsl)
inputf.close()

mol, atids, resids = get_atominfo_fpdb(options.pdbf)
mol = charmm_mol(mol, atids)

if len(atids) != max(atids):
    raise ValueError('The atom numbers are not consecutive!')

for i in bpairsl:
    w_output = open(i[2], 'w')
    for j in xrange(i[0], i[1]+1):
        for k in mol.residues[j].resconter:
            gtype = mol.atoms[k].gtype
            atid = mol.atoms[k].atid
            atname = mol.atoms[k].atname
            crd = mol.atoms[k].crd
            resid = mol.atoms[k].resid
            resname = mol.atoms[k].resname
            if len(atname) == 3:
                print("%-6s%5d %4s %3s %1s%4d   %8.3f%8.3f%8.3f%6.2f%6.2f" %(gtype, atid, atname,
                resname, ' ', resid, crd[0], crd[1], crd[2], 1.00, 0.00), file=w_output)
            else:
                print("%-6s%5d %4s %3s %1s%4d   %8.3f%8.3f%8.3f%6.2f%6.2f" %(gtype, atid, atname.center(4),
                resname, ' ', resid, crd[0], crd[1], crd[2], 1.00, 0.00), file=w_output)
    w_output.close()
    if options.renum == 1:
        os.system("mv %s %s" %(i[2], i[2]+'.bak'))
        os.system("pdb4amber -i %s -o %s" %(i[2]+'.bak', i[2]))
    fname = i[2]
    charmmf = fname.rstrip('.pdb')
    charmmf = charmmf + '.in'
    write_charmm_inp(fname, charmmf, options.rtf, options.prm)

if options.water != 0:
    waterl = []
    for i in resids:
        resname = mol.residues[i].resname
        if resname in ['WAT', 'HOH', 'TP3', 'TIP3']:
            waterl.append(i)
    for i in bpairsl:
        if min(waterl) <= i[1]:
            raise ValueError('Water residues should not contained in '
                'the input file ranges if -w equals 1!')    
    nwater1 = len(waterl)
    nwater2 = max(waterl)-min(waterl)+1
    if nwater1 != nwater2:
        raise ValueError('Water residues are not consecutive!')
    else:
        # Cut the water residues into different residues
        if nwater1%options.water == 0:
            nfiles = nwater1/options.water
        else:
            nfiles = nwater1/options.water + 1

        for i in xrange(1, nfiles+1):
            tmpl = []
            fname = 'water' + str(i) + '.pdb'
            w_output = open(fname, 'w')
            for j in range((i-1)*options.water+min(waterl), i*options.water+min(waterl)):
                if (j in resids):
                    for k in mol.residues[j].resconter:
                        gtype = mol.atoms[k].gtype
                        atid = mol.atoms[k].atid
                        atname = mol.atoms[k].atname
                        crd = mol.atoms[k].crd
                        resid = mol.atoms[k].resid
                        resname = mol.atoms[k].resname
                        if len(resname) == 3 and len(atname) == 3:
                            print("%-6s%5d %4s %3s %1s%4d   %8.3f%8.3f%8.3f%6.2f%6.2f" %(gtype, atid, atname,
                            resname, ' ', resid, crd[0], crd[1], crd[2], 1.00, 0.00), file=w_output)
                        elif len(resname) == 3 and len(atname) in [2, 4]:
                            print("%-6s%5d %4s %3s %1s%4d   %8.3f%8.3f%8.3f%6.2f%6.2f" %(gtype, atid, atname.center(4),
                            resname, ' ', resid, crd[0], crd[1], crd[2], 1.00, 0.00), file=w_output)
                        elif len(resname) == 4 and len(atname) == 3:
                            print("%-6s%5d %4s %4s%1s%4d   %8.3f%8.3f%8.3f%6.2f%6.2f" %(gtype, atid, atname,
                            resname, ' ', resid, crd[0], crd[1], crd[2], 1.00, 0.00), file=w_output)
                        elif len(resname) == 4 and len(atname) in [2, 4]:
                            print("%-6s%5d %4s %4s%1s%4d   %8.3f%8.3f%8.3f%6.2f%6.2f" %(gtype, atid, atname.center(4),
                            resname, ' ', resid, crd[0], crd[1], crd[2], 1.00, 0.00), file=w_output)
            tmpl.append(j)
            w_output.close()
            bpairsl.append([min(tmpl), max(tmpl), fname])
            if options.renum == 1:
                bakname = fname + '.bak'
                os.system("mv %s %s" %(fname, bakname))
                os.system("pdb4amber -i %s -o %s" %(bakname, fname))
            charmmf = 'water' + str(i) + '.in'
            write_charmm_inp(fname, charmmf, options.rtf, options.prm)

# Write the CHARMM input file
write_charmm_combf(bpairsl, 'combined_react.in', options.rtf, options.prm)
write_charmm_combf(bpairsl, 'combined_product.in', options.rtf, options.rtf)

