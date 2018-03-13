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

import Bio
import Bio.PDB
import numpy
from numpy import *

from Bio.SVDSuperimposer import SVDSuperimposer

pdbs= []
inp=open(sys.argv[3], "r")

while 1:
	line = inp.readline()
    	if not line: break
	line=line.split()
	pdbs.append((str(line[0]),str(line[2])))

print pdbs

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

#given a fixed and moving set of atoms, this calculates the local rmsd for each residue
def localrmsd(fixed, moving):
	lrms = []
	if not (len(fixed)==len(moving)):
	    raise PDBException("Fixed and moving atom lists differ in size")
	l=len(fixed)
	fixed_coord=numpy.zeros((l, 3))
	moving_coord=numpy.zeros((l, 3))
	for i in range(0, len(fixed)):
	    fixed_coord[i]=fixed[i].get_coord()
	    moving_coord[i]=moving[i].get_coord()
	    diff=fixed_coord[i] - moving_coord[i]
	    lrms.append(sqrt(sum(diff*diff)))
	return lrms

#method to calculate the local rmsds for the given range0 to range1 from the structure stru and then save it as bfactors in model mod
def conf(stru,listofres, mod):
	ave=[]
	nmod=len(stru)
	for i, refm in zip(range(0,nmod), stru):
		for j, altm in zip(range(0,nmod), stru) :
		#Build paired lists of c-alpha atoms:
			rloop = []
			aloop = []
			for (refc, altc) in zip(refm, altm) :
				for refr, altr in zip(refc, altc) :
					assert refr.resname == altr.resname
					assert refr.id      == altr.id
					#if refr.id[1] > range0 and refr.id[1] < range1:
					if refr.id[1] in listofres:
					    rloop.append(refr['CA'])                
					    aloop.append(altr['CA'])

			ave.append(localrmsd(rloop, aloop))
	stdevs=numpy.std(array(ave),axis=0)
	for c in mod:
		for r in c:
			if r.id[1] in listofres:
				for a in r:
					a.set_bfactor(stdevs[listofres.index(r.id[1])])
			else:
				for a in r:
					a.set_bfactor(0.0)
	
	return stdevs


#pymol string
pymolstr='load\ \%s,prot\;hide\ all\;show\ cartoon\;spectrum\ b\;cartoon\ putty\;viewport\ 300,300\;ray\ 800,800\;png\ \%s\.png\;orient\;mset\ 1\ x60\;mplay\;util.mroll\ 1,60,180\;mstop\;set\ ray_trace_frames=1\;mpng\ frame'


#pdb file
pdb_code = sys.argv[1]
pdb_filename = "%s.pdb" % pdb_code


print "Loading PDB file %s" % pdb_filename
structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)

nmodels=len(structure)
mat=numpy.zeros((nmodels,nmodels))

loc=[]

inp2=open(sys.argv[2], "r")
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

#This aligns all the structures (all but loop)
ref_model=structure[0]
print "Everything aligned to first model.."
for j, alt_model in zip(range(0,nmodels), structure) :
   ref_atoms = []
   alt_atoms = []
   for (ref_chain, alt_chain) in zip(ref_model, alt_model) :
	   for ref_res, alt_res in zip(ref_chain, alt_chain) :
		assert ref_res.resname == alt_res.resname
		assert ref_res.id      == alt_res.id
		#if ref_res.id[1] < mini or ref_res.id[1] > maxi:
		if not ref_res.id[1] in loopres:
		    ref_atoms.append(ref_res['CA'])                
		    alt_atoms.append(alt_res['CA'])

		#Align these paired atom lists:
		super_imposer = Bio.PDB.Superimposer()
		super_imposer.set_atoms(ref_atoms, alt_atoms)

	   if ref_model.id == alt_model.id :
		assert numpy.abs(super_imposer.rms) < 0.0000001
	    	assert numpy.max(numpy.abs(super_imposer.rotran[1])) < 0.000001
	    	assert numpy.max(numpy.abs(super_imposer.rotran[0]) - numpy.identity(3)) < 0.000001
    	   else :
	   	 super_imposer.apply(alt_model.get_atoms())

#After the alignment, calculate the loop rmsd
for i, ref_model in zip(range(0,nmodels), structure):
	for j, alt_model in zip(range(0,nmodels), structure) :
		#Build paired lists of c-alpha atoms:
		ref_loop = []
		alt_loop = []
		for (ref_chain, alt_chain) in zip(ref_model, alt_model) :
			for ref_res, alt_res in zip(ref_chain, alt_chain) :
				assert ref_res.resname == alt_res.resname
				assert ref_res.id      == alt_res.id
				#if ref_res.id[1] > mini and ref_res.id[1] < maxi:
				if ref_res.id[1] in loopres:
			            ref_loop.append(ref_res['CA'])                
				    alt_loop.append(alt_res['CA'])

		mat[i][j] = round(set_atoms2(ref_loop,alt_loop),2)

#print "Saving aligned structure as PDB file %s" % pdb_out_filename
io=Bio.PDB.PDBIO()
#io.set_structure(structure)
#io.save(pdb_out_filename)

print mat

#writing rmsd matrix to file
outfile=open('rmsdout', 'w')
for i in range(nmodels): 
	for j in range(nmodels):
		outfile.write(str(mat[i][j]) + '\t')
        outfile.write('\n')
outfile.close()

print "Done generating matrix"


#import Numeric
from Bio import Cluster

matt=array(mat)

#Clustering
tree= Cluster.treecluster(data=None,method='m',  distancematrix=matt)

print "Tree is ", tree

#change this cutat value to set the cutoff for the tree
cutat=float(sys.argv[4])
count=0
'''
for i in range(1,10):
	for node in tree:
		if(node.distance > i):
			count = count + 1
	print "The number of clusters for %s cutoff is %s" %(i,count)
	count=0
'''

#finding the number of clusters for the given cutat
for node in tree:
	if(node.distance > cutat):
		count = count + 1

print "The total number of cluster for cutoff of %s A is %s " %(cutat,count) 

#cid stores the cluster id for each element 
cid=tree.cut(count)
print "cid is ", cid

#arrc stores the info about cluster and its size
arrc=[]
for i in range(0,count):
	print "Size of cluster %s is %s" %(i,list(cid).count(i))
	arrc.append((list(cid).count(i),i))

#largest cluster info is in arrc[0] 
print "arr is ",arrc
arrc.sort(reverse=True)	
print "arr sorted is ", arrc

l=0
bigg=0
topc=[[] for ni in range(count)]

#topc stores cluster info. For eg. topc[0] contains an array of ids of those elements that belong to the top cluster (Cluster 0)		
for i in range(count):
	for j,item in enumerate(cid):
		if(item==arrc[i][1]):
			topc[i].append(j)
#print "The biggest cluster is %s with size %s" %(bigg,l)
#print topc.count()
print topc


cent=numpy.zeros(len(topc))
#finding the centroids of each cluster (centroid is simply defined here as the elemet in the cluster whose total distance to all other elements in that cluster is the smallest
#cent2ave and cent2std contain the average rmsd of the cluster and the average stdev for that cluster
cent2=numpy.zeros(len(topc))
cent2ave=numpy.zeros(len(topc))
cent2std=numpy.zeros(len(topc))
clusaveE=numpy.zeros(len(topc))
for i in range(0,len(topc)):
	temp=-1
	topcc=numpy.zeros((len(topc[i]),len(topc[i])))
	aveE=numpy.zeros(len(topc[i]))
	for j,ji in enumerate(topc[i]):
		for k,ki in enumerate(topc[i]):
			topcc[j][k]=matt[ji][ki]
		if(cent[i] > numpy.sum(topcc[j]) or temp < 0 ):
				cent[i] = numpy.sum(topcc[j])
				cent2[i] = ji
				temp=cent[i]
		print ji, numpy.sum(topcc[j])
		aveE[j]=pdbs[ji][1]

	cent2ave[i]=numpy.average(topcc)
	cent2std[i]=numpy.std(topcc)
	clusaveE[i]=numpy.average(aveE)
print "Centroid is ", cent2


#cluster data file
clusdata=open('clusterdata.csv', 'w')

clusdata.write("Cluster,Ave RMSD, Stdev RMSD,Ave Energy,Cluster Size\n");

#saving the top 5 clusters along with the centroids with bfactor as local rmsd and output a pymol movie as well with the local rmsd as confidence
ntopclus=5
#a=[]
from Bio.PDB.Structure import Structure
for i in range(ntopclus):
	s=Structure(0)
	c=Structure(0)
	for j in topc[i]:
		print j 
		s.add(structure[j])
	m=structure[cent2[i]]
	print m
	print conf(s,loopres,m)
	c.add(m)
	pdb_out = pdb_code +  "Cluster" + str(i) + ".pdb" 
	cent_out = pdb_code +  "Cluster" + str(i) + "Centroid.pdb" 
	print "The centroid of cluster %s is %s and the ave rmsd of the cluster is %s with a stdev of %s and average energy of %s " %(str(i), str(cent2[i]), str(cent2ave[i]),str(cent2std[i]),str(clusaveE[i]) )
	print "Saving aligned structure as PDB file %s" % pdb_out
	io.set_structure(s)
	io.save(pdb_out)
	io.set_structure(c)
	io.save(cent_out)
	os.system('pymol -c -d ' + (pymolstr %(cent_out, cent_out)))
	os.system('convert frame*.png %s.gif' %(cent_out))
	clusdata.write(str(i) + ',' + str(round(cent2ave[i],2)) +  ',' + str(cent2std[i]) + ','  + str(clusaveE[i]) + ',' + str(len(topc[i])) + '\n')

clusdata.close()

#Make gnuplot plot for cluster data
'''
set term png\n
set samples 5000\n
set output "ClusterStat%s.png"\n
set title "Cluster Statistics"\n
set ylabel "Average RMSD within cluster (A)"\n
set xlabel "Cluster Number"\n
set cblabel "Energy"\n
set xrange [-1:5]\n
set view map\n
set hidden3d\n
splot 'clusterdata' u 1:2:4:5 w p pt 7 ps var lt palette notitle, '' u 1:2:(0):5 w labels notitle
'''


#when native known
'''
set term png\n
set samples 5000\n
set output "ClusterStat.png"\n
set title "Cluster Statistics for the five largest loop clusters for target T0411D1"\n
set ylabel "Cluster Size"\n
set xlabel "DOPE-PW Energy"\n
set cblabel "Tightness (Average RMSD within cluster)"\n
set xrange [-2000:-1600]\n
set yrange [5:10]\n
set view map\n
#set hidden3d\n
set palette rgbformulae 32,13,10
#set palette rgbformulae 23,28,3
splot 'clusterdata' u 4:5:2:($6) w p pt 7 ps var lt palette notitle,'' u 4:5:(0):6 w labels notitle\n
'''

gnuplotstr='''
set term png\n
set samples 5000\n
set output "ClusterStat.png"\n
set title "Cluster Statistics for the five largest loop clusters"\n
set ylabel "Cluster Size"\n
set xlabel "DOPE-PW Energy"\n
set cblabel "Tightness (Average RMSD within cluster)"\n
#set xrange [-2000:-1600]\n
#set yrange [5:10]\n
set view map\n
#set hidden3d\n
set palette rgbformulae 32,13,10
#set palette rgbformulae 23,28,3
splot '<sed "s/,/\ /g" clusterdata' u 4:5:2:($2) w p pt 7 ps var lt palette notitle,'' u 4:5:(0):2 w labels notitle\n
'''
 
print "Gnuplotting .... "

gnuscript=open('gnufile', 'w')
#gnuscript.write(gnuplotstr % (pdb_code))
gnuscript.write(gnuplotstr)
#this is gnuplot 4.3 installed in /home/aashish/gnuplot in CI
gnuscript.close()
os.system('/home/aashish/gnuplot/bin/gnuplot gnufile')
