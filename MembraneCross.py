#!/usr/bin/env python

import sys

from Bio.PDB import *
from Bio.PDB.Fold import *

#structure is arg1, everything must be chain A
stru=PDBParser().get_structure('x',sys.argv[1])


#A,B,C define plane (CA residues). D is the top of the loop. E is inside the
#membrane
A=int(sys.argv[2])
B=int(sys.argv[3])
C=int(sys.argv[4])
D=int(sys.argv[5])
E=int(sys.argv[6])


import numpy

a=numpy.array(stru[0]['A'][A]['CA'].get_coord())
b=numpy.array(stru[0]['A'][B]['CA'].get_coord())
c=numpy.array(stru[0]['A'][C]['CA'].get_coord())
d=numpy.array(stru[0]['A'][D]['CA'].get_coord())
e=numpy.array(stru[0]['A'][E]['CA'].get_coord())

bp=b-a
cp=c-a
dp=d-a
ep=e-a

side1=numpy.linalg.det(numpy.array([bp,cp,dp]))
side2=numpy.linalg.det(numpy.array([bp,cp,ep]))


def sign(number):return cmp(number,0)

#if 1, same side, if 0 different side
if sign(side1)==sign(side2):
    print 1
else:
    print 0

