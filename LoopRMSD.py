#Clustering script for loop modeling results. works for multiple ranges
#.A.

#Usage: python script.py arg1 arg2 arg3 arg4
# arg1 is the pdb code xxxx (should have a file called xxxx.pdb in the folder that contains all the final pdb structures) 
# arg2 is the range file (relative to first res)
# arg3 is the file containing the list of pdbs and energies
# arg4 is the cutoff value to cut the tree
# 
# arg3 is of the form:
#			col1	 col2	  col3
#			xxxx.pdb ligation dopepw
#
# Note: The clustering in this is based on global loop rmsd - ie. align nonloop and calculate rmsd of loop residues

import sys,os
sys.path.insert(0, '/home/aashish/Loops/lib-loops/python2.6/site-packages')

import Bio
import Bio.PDB
import numpy
from numpy import *

from Bio.SVDSuperimposer import SVDSuperimposer

#this calculates the rmsd between the given set of fixed and moving - (sum of the local rmsds)
def set_atoms2(fixed, moving):
	"""
	Put (translate/rotate) the atoms in fixed on the atoms in 
	moving, in such a way that the RMSD is minimized.

	@param fixed: list of (fixed) atoms
	@param moving: list of (moving) atoms 
	@type fixed,moving: [L{Atom}, L{Atom},...]
	"""
	if not (len(fixed)==len(moving)):
		 raise PDBException("Fixed and moving atom lists differ in size")
	l=len(fixed)
	fixed_coord=numpy.zeros((l, 3))
	moving_coord=numpy.zeros((l, 3))
	for i in range(0, len(fixed)):
		fixed_coord[i]=fixed[i].get_coord()
		moving_coord[i]=moving[i].get_coord()

	diff=fixed_coord-moving_coord
	l=fixed_coord.shape[0]
	return sqrt(sum(sum(diff*diff))/l)


#pdb file
native_file = sys.argv[1]
pdt_file = sys.argv[2]

inp2=open(sys.argv[3], "r")
loopres=[]
while 1:
	line = inp2.readline()
    	if not line: break
	if "NR" in line: 
		break
	line=line.split()
	mini=int(line[0])
	maxi=int(line[1])
	for i in range(mini,maxi+1):
		loopres.append(i)

print " Loop Residues are: ", loopres

dolocal=0
doglobal=0
dowhat=sys.argv[4]
if dowhat=="local": dolocal=1
if dowhat=="global": doglobal=1

print dolocal, doglobal

print "Loading PDB file %s" % native_file
structure = Bio.PDB.PDBParser().get_structure('native', native_file)
structure2 = Bio.PDB.PDBParser().get_structure('others', pdt_file)


ref_model=structure[0]

nmodels=len(structure2)
#mat=numpy.zeros((nmodels,nmodels))


#This aligns all the structures (all but loop)
print "Everything aligned to the native model.."
for j, alt_model in zip(range(0,nmodels), structure2) :
   ref_atoms = []
   alt_atoms = []
   ref_loop = []
   alt_loop = []
	
   for (ref_chain, alt_chain) in zip(ref_model, alt_model) :
	   for ref_res, alt_res in zip(ref_chain, alt_chain) :
		assert ref_res.resname == alt_res.resname
		assert ref_res.id      == alt_res.id
		#if ref_res.id[1] < mini or ref_res.id[1] > maxi:
		
		#this makes sure the rmsd is for the loop. Remove this condition for global rmsd between two structures
		if ref_res.id[1] in loopres:
		    ref_loop.append(ref_res['CA'])                
		    alt_loop.append(alt_res['CA'])
		

		#these make sure what to align
		if not (ref_res.id[1] in loopres) and doglobal:
		    ref_atoms.append(ref_res['CA'])                
		    alt_atoms.append(alt_res['CA'])
		if ref_res.id[1] in loopres and dolocal:
		    ref_atoms.append(ref_res['CA'])                
		    alt_atoms.append(alt_res['CA'])
#		else:
#		    ref_atoms.append(ref_res['CA'])                
#		    alt_atoms.append(alt_res['CA'])

	#Align these paired atom lists:
   super_imposer = Bio.PDB.Superimposer()
   super_imposer.set_atoms(ref_atoms, alt_atoms)

#	   if ref_model.id == alt_model.id :
#		 print "Ignoring self self RMS"
	#assert numpy.abs(super_imposer.rms) < 0.0000001
	#assert numpy.max(numpy.abs(super_imposer.rotran[1])) < 0.000001
	#assert numpy.max(numpy.abs(super_imposer.rotran[0]) - numpy.identity(3)) < 0.000001
#   	   else :
   super_imposer.apply(alt_model.get_atoms())
   print "Native and Model %i - Align: %0.2f  RMSD: %0.2f" % (alt_model.id, super_imposer.rms, set_atoms2(ref_loop,alt_loop))

pdb_out_filename = "%s_aligned.pdb" % pdt_file
io=Bio.PDB.PDBIO()
io.set_structure(structure2)
io.save(pdb_out_filename)

print "Saving aligned structure as PDB file %s" % pdb_out_filename
