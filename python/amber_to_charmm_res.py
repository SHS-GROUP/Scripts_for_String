#!/home/pengfeil/AMBER/amber16/amber16/miniconda/bin/python
# Filename: amber_to_charmm_res.py

###############################################################################
#                            FOR NUCLEIC ACIDS
###############################################################################


NA_ATOM_DICT = {'H5\'1' : 'H5\'',
                'H5\'2' : 'H5\'\'',
                'H2\'1' : 'H2\'\'',
                'HO\'2' : 'H2\'',
                '\'HO3' : 'H3T',
                'HO5\'' : 'H5T',
                'OXT'   : 'OT2',
               }

NA_RES_DICT = {
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

###############################################################################
#                            FOR AMINO ACIDS
###############################################################################

# In total there are 28 normal AAs, 24 N-temrinal AAs (containing ACE),
# and 26 C-terminal AAs (containg NHE and NME) in the AMBER ff14SB
# for protein

AA_RES_LIST = {'ALA', 'ARG', 'ASH', 'ASN', 'ASP', #There is no NASH, CASH
               'CYM', 'CYS', 'CYX',  #There is no NCYM, CCYM
               'GLH', 'GLN', 'GLU', 'GLY', #There is no NGLH, CGLH
               'HID', 'HIE', 'HIP', 'HYP', #There is no NHYP
               'ILE', 'LEU', 'LYN', 'LYS', #There is no NLYN, CLYN
               'MET', 'PHE', 'PRO', 'SER',
               'THR', 'TRP', 'TYR', 'VAL',
               'ACE', 'NHE', 'NME'}

# ASH is protonated ASP
# GLH is protonated GLU
# CYM is unprotonated CYS
# CYX is S-S bonded CYS
# HID, HIE, and HIP are HISs with different protonated states
# HYP is Hydro-PRO
# LYN is unprotonated LYS
# ACE is CH3-CO
# NHE is NH2
# NME is CH3-NH

AA_RES_ATOM_DICT = {
	('ALA', 'H') : ('ALA', 'HN'),     #ALA
	('ARG', 'H') : ('ARG', 'HN'),     #ARG
	('ARG', 'HB2') : ('ARG', 'HB1'), 
        ('ARG', 'HB3') : ('ARG', 'HB2'),
	('ARG', 'HG2') : ('ARG', 'HG1'),
        ('ARG', 'HG3') : ('ARG', 'HG2'),
        ('ARG', 'HD2') : ('ARG', 'HD1'),
        ('ARG', 'HD3') : ('ARG', 'HD2'),
        # There is no ASH in CHARMM parm14sb_all.rtf
     	('ASN', 'H')  : ('ASN', 'HN'),    #ASN
	('ASN', 'HB2') : ('ASN', 'HB1'),
        ('ASN', 'HB3') : ('ASN', 'HB2'),
	('ASP', 'H') : ('ASP', 'HN'),     #ASP
	('ASP', 'HB2') : ('ASP', 'HB1'),
        ('ASP', 'HB3') : ('ASP', 'HB2'),
	('CYM', 'H') : ('CYM', 'HN'),     #CYM
	('CYM', 'HB2') : ('CYM', 'HB1'),
	('CYM', 'HB3') : ('CYM', 'HB2'),
	('CYS', 'H') : ('CYS', 'HN'),     #CYS
        ('CYS', 'HB2') : ('CYS', 'HB1'),
        ('CYS', 'HB3') : ('CYS', 'HB2'),
	('CYS', 'HG') : ('CYS', 'HG1'),
	('CYX', 'H') : ('CYX', 'HN'),     #CYX
        ('CYX', 'HB2') : ('CYX', 'HB1'),
        ('CYX', 'HB3') : ('CYX', 'HB2'),
	('GLH', 'H') : ('GLH', 'HN'),     #GLH
	('GLH', 'HB2') : ('GLH', 'HB1'),
	('GLH', 'HB3') : ('GLH', 'HB2'),
    	('GLH', 'HG2') : ('GLH', 'HG1'),
    	('GLH', 'HG3') : ('GLH', 'HG2'),
	('GLN', 'H') : ('GLN', 'HN'),     #GLN
	('GLN', 'HB2') : ('GLN', 'HB1'),
	('GLN', 'HB3') : ('GLN', 'HB2'),
    	('GLN', 'HG2') : ('GLN', 'HG1'),
    	('GLN', 'HG3') : ('GLN', 'HG2'),
    	('GLU', 'H') : ('GLU', 'HN'),     #GLU
        ('GLU', 'HB2') : ('GLU', 'HB1'),
        ('GLU', 'HB3') : ('GLU', 'HB2'),
        ('GLU', 'HG2') : ('GLU', 'HG1'),
        ('GLU', 'HG3') : ('GLU', 'HG2'),
    	('GLY', 'H') : ('GLY', 'HN'),     #GLY
    	('GLY', 'HA2') : ('GLY', 'HA1'),
        ('GLY', 'HA3') : ('GLY', 'HA2'),
    	('HID', 'H') : ('HID', 'HN'),     #HID
    	('HID', 'HB2') : ('HID', 'HB1'),
        ('HID', 'HB3') : ('HID', 'HB2'),
    	('HIE', 'H') : ('HIE', 'HN'),     #HIE
    	('HIE', 'HB2') : ('HIE', 'HB1'),
        ('HIE', 'HB3') : ('HIE', 'HB2'),
    	('HIP', 'H') : ('HIP', 'HN'),     #HIP
    	('HIP', 'HB2') : ('HIP', 'HB1'),
        ('HIP', 'HB3') : ('HIP', 'HB2'),
        # There is no HYP in CHARMM parm14sb_all.rtf
    	('ILE', 'H') : ('ILE', 'HN'),     #ILE
    	('ILE', 'CD1') : ('ILE', 'CD'),
    	('ILE', 'HD11') : ('ILE', 'HD1'),
    	('ILE', 'HD12') : ('ILE', 'HD2'),
    	('ILE', 'HD13') : ('ILE', 'HD3'),
    	('ILE', 'HG12') : ('ILE', 'HG11'),
        ('ILE', 'HG13') : ('ILE', 'HG12'),
        ('LEU', 'H') : ('LEU', 'HN'),     #LEU
        ('LEU', 'HB2') : ('LEU', 'HB1'),
        ('LEU', 'HB3') : ('LEU', 'HB2'),
        # There is no LYN in CHARMM parm14sb_all.rtf
    	('LYS', 'H') : ('LYS', 'HN'),     #LYS
    	('LYS', 'HB2') : ('LYS', 'HB1'),
        ('LYS', 'HB3') : ('LYS', 'HB2'),
    	('LYS', 'HG2') : ('LYS', 'HG1'),
        ('LYS', 'HG3') : ('LYS', 'HG2'),
    	('LYS', 'HD2') : ('LYS', 'HD1'),
        ('LYS', 'HD3') : ('LYS', 'HD2'),
    	('LYS', 'HE2') : ('LYS', 'HE1'),
        ('LYS', 'HE3') : ('LYS', 'HE2'),
    	('MET', 'H') : ('MET', 'HN'),     #MET
        ('MET', 'HB2') : ('MET', 'HB1'),
        ('MET', 'HB3') : ('MET', 'HB2'),
    	('MET', 'HG2') : ('MET', 'HG1'),
        ('MET', 'HG3') : ('MET', 'HG2'),
    	('PHE', 'H') : ('PHE', 'HN'),     #PHE
        ('PHE', 'HB2') : ('PHE', 'HB1'),
        ('PHE', 'HB3') : ('PHE', 'HB2'),
    	('PRO', 'HD2') : ('PRO', 'HD1'),  #PRO
        ('PRO', 'HD3') : ('PRO', 'HD2'),
    	('PRO', 'HG2') : ('PRO', 'HG1'),
        ('PRO', 'HG3') : ('PRO', 'HG2'),
    	('PRO', 'HB2') : ('PRO', 'HB1'),
        ('PRO', 'HB3') : ('PRO', 'HB2'),
    	('SER', 'H') : ('SER', 'HN'),     #SER
    	('SER', 'HB2') : ('SER', 'HB1'),
        ('SER', 'HB3') : ('SER', 'HB2'),
    	('SER', 'HG') : ('SER', 'HG1'),
    	('THR', 'H') : ('THR', 'HN'),     #THR
    	('TRP', 'H') : ('THR', 'HN'),     #TRP
    	('TRP', 'HB2') : ('TRP', 'HB1'),  
        ('TRP', 'HB3') : ('TRP', 'HB2'),
    	('TYR', 'H') : ('TYR', 'HN'),     #TYR
    	('TYR', 'HB2') : ('TYR', 'HB1'),
        ('TYR', 'HB3') : ('TYR', 'HB2'),
    	('VAL', 'H') : ('VAL', 'HN'),     #VAL
        # There is no NHE in CHARMM parm14sb_all.rtf
        ('NME', 'H') : ('NME', 'HN'),     #NME
        ('NME', 'CH3') : ('NME', 'CAT'),
        ('NME', 'HH31') : ('NME', 'HT1'),
        ('NME', 'HH32') : ('NME', 'HT2'),
        ('NME', 'HH33') : ('NME', 'HT3'),
        ('ACE', 'CH3') : ('ACE', 'CAY'),  #ACE
        ('ACE', 'HH31') : ('ACE', 'HY1'),
        ('ACE', 'HH32') : ('ACE', 'HY2'),
        ('ACE', 'HH33') : ('ACE', 'HY3'),
        ('WAT', 'O')   : ('TIP3', 'OH2'), #WAT
        ('WAT', 'H1')  : ('TIP3', 'H1'),
        ('WAT', 'H2')  : ('TIP3', 'H2'),  
    	('Na+', 'Na+') : ('SOD', 'SOD'),  #Na+
    	('Cl-', 'Cl-') : ('CLA', 'CLA'),  #Cl-
}

TER_AA_RES_ATOM_DICT = {} #Terminal AA atom dict
for aa in AA_RES_LIST:
    if aa not in ['ASH', 'CYM', 'GLH', 'HYP', 'LYN', 'ACE', 'NHE', 'NME']:
        TER_AA_RES_ATOM_DICT[(aa,'H1')] = ('N' + aa, 'HT1')
        TER_AA_RES_ATOM_DICT[(aa,'H2')] = ('N' + aa, 'HT2')
        TER_AA_RES_ATOM_DICT[(aa,'H3')] = ('N' + aa, 'HT3')
    if aa not in ['ASH', 'CYM', 'GLH', 'LYN', 'ACE', 'NHE', 'NME']:
        TER_AA_RES_ATOM_DICT[(aa,'O')] = ('C' + aa, 'OT1')
        TER_AA_RES_ATOM_DICT[(aa,'OXT')] = ('C' + aa, 'OT2')

###############################################################################
#                            FOR ATOM TYPES IN FF14SB
###############################################################################

"""
ATOM_TYPE_DICT = {'c' : 'C',
                  'c3': 'CT',
                  'o' : 'O2',
                  'c2': 'CM',
                  'hc': 'HC',
                  'ha': 'HA'}
"""

ATOM_TYPE_DICT = {'Na+': 'SOD',
                  'Cl-': 'CLA',
                  'C*' : 'CG',
                  'N*' : 'NG'
                 }


