#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
import parmed as pmd
from optparse import OptionParser
import os
from amber_to_charmm_res import (AA_RES_LIST, AA_RES_ATOM_DICT,
                          TER_AA_RES_ATOM_DICT, ATOM_TYPE_DICT)

parser = OptionParser("Usage: amb2chm_psf_crd.py -p prmtop -c inpcrd -f psf -d crd -b pdb")
parser.add_option("-p", dest="prmtop", type='string',
                  help="Prmtop file")
parser.add_option("-c", dest="inpcrd", type='string',
                  help="Inpcrd file")
parser.add_option("-f", dest="psf", type='string',
                  help="Psf file")
parser.add_option("-d", dest="crd", type='string',
                  help="Crd file")
parser.add_option("-b", dest="pdb", type='string',
                  help="Crd file")
(options, args) = parser.parse_args()

amber = pmd.load_file(options.prmtop, options.inpcrd)

# For the normal residues and atom types
for i in xrange(0, len(amber.atoms)):
    resname = amber.atoms[i].residue.name
    atname = amber.atoms[i].name
    attype = amber.atoms[i].type
    if (resname, atname) in AA_RES_ATOM_DICT:
        amber.atoms[i].residue.name = AA_RES_ATOM_DICT[(resname, atname)][0]
        amber.atoms[i].name = AA_RES_ATOM_DICT[(resname, atname)][1]
    if attype in ATOM_TYPE_DICT:
        amber.atoms[i].type = ATOM_TYPE_DICT[attype]

# For the N-temrinal and C-terminal amino acids
for i in xrange(0, len(amber.residues)):
    resname = amber.residues[i].name
    if resname in AA_RES_LIST:
        atnames = []
        for j in xrange(0, len(amber.residues[i].atoms)):
            atname = amber.residues[i].atoms[j].name
            atnames.append(atname)
        if 'H1' in atnames and 'H2' in atnames and 'H3' in atnames:
            amber.residues[i].name = 'N' + resname
            for j in xrange(0, len(amber.residues[i].atoms)):
                atname = amber.residues[i].atoms[j].name
                if atname == 'H1':
                    amber.residues[i].atoms[j].name = 'HT1'
                elif atname == 'H2':
                    amber.residues[i].atoms[j].name = 'HT2'
                elif atname == 'H3':
                    amber.residues[i].atoms[j].name = 'HT3'
        elif 'O' in atnames and 'OXT' in atnames:
            amber.residues[i].name = 'C' + resname
            for j in xrange(0, len(amber.residues[i].atoms)):
               atname = amber.residues[i].atoms[j].name
               if atname == 'O':
                   amber.residues[i].atoms[j].name = 'OT1'
               elif atname == 'OXT':
                   amber.residues[i].atoms[j].name = 'OT2'

# Save a CHARMM PSF and crd file
amber.save(options.psf, overwrite=True, format='psf')
amber.save(options.crd, overwrite=True, format='charmmcrd')
amber.save(options.pdb, overwrite=True, format='pdb')

quit()

