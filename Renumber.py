#!/usr/bin/env python

import sys

from Bio.PDB import *
from Bio.PDB.Fold import *

stru=PDBParser().get_structure('x',sys.argv[1])

i=int(sys.argv[2])

for model in stru:
    for chain in model:
        for res in chain:
            res.id = (' ',i, ' ')
            i += 1
w = PDBIO()
w.set_structure(stru)
renum=sys.argv[1] + '-renum.pdb'
w.save(renum)
