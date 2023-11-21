#!/usr/bin/python

#Weihua Zheng on Mar,2017
#Take a PDB file as an input, calculate (1-Q_relative) for neighboring snapshots
#This code deals with specifically the case when each chain is the same as others
#The scrambling of the chains creates a permutation problem, such as in the case
#when two essentially identical structures have different chain orders, they will
#have very low Q-relative, instead of a high Q-relative.  
# ----------------------------------------------------------------------

import sys
#from skimage._build import cython
if len(sys.argv)!=7 :	
	print "\n exe PDBs-Input Output_Q-relative n_chains permutation_flag(0/1) Qonuchic_flag(0/1) n_processor(1/2/3/4) \n"
	print "***************************"
	print "Q-rel pair is in the following order: 1 2, 1 3, ...,1 n; 2 3, 2 4, ..., 2 n; ...; n-1 n."
	print "It is a vector with size n*(n-1)/2"
	print "Uses Qcontact when Q_onuchic flag is on "
	sys.exit()

import Zheng_func
from itertools import permutations

#A modified version of ComputeQ_permutation() for multiprocessor
#return 1-Q instead of Q
def computeD_multi_processor(lists_cb_atoms):
	cb_atoms_1 = lists_cb_atoms[0]
	cb_atoms_2 = lists_cb_atoms[1]
	N  = len(cb_atoms_1)
	N2 = len(cb_atoms_2)
	if N != N2:
		print "Error, cb_atoms length mismatch."
		print "Length of cb_atoms_1 is ", N, "Length of cb_atoms_2 is ", N2
		print "Check the number of residues in the pdb files."
		sys.exit()
	
	if perm_flag == 0: #global variable
		Q_rel = Zheng_func.computeQ_relative(cb_atoms_1, cb_atoms_2, qo_flag, cutoff)

	#create all permuated cb_atoms_2_perm list for Q_rel calculation, max Q_rel is selected.
	else:		
		Q_perm = []
		n_res = N / n_chains #number of res per chain
		assert(n_res == int(N/n_chains)),"Each chain should have the same # of residues."
		
		#cool stuff		
		residue_reshape_list = [range(N)[i:i+n_res] for i in range(0,N,n_res)]
		
		residue_list = range(N)
		chain_list   = range(n_chains)					
		n_chain_perm = list(permutations(chain_list))
		
		for chainlist_i in n_chain_perm:
			residue_list_perm = []
			for id_chain in chainlist_i :				
				residue_list_perm.extend(residue_reshape_list[id_chain])

			cb_atoms_2_perm = [cb_atoms_2[residue_list_perm[id_res]] for id_res in residue_list]

			Q = Zheng_func.computeQ_relative(cb_atoms_1, cb_atoms_2_perm, qo_flag, cutoff)
			Q_perm.append(Q)
		Q_rel = max(Q_perm)
	return 1 - Q_rel

pdb_file    = sys.argv[1]
output_file = sys.argv[2]
n_chains    = int(sys.argv[3])
perm_flag   = int(sys.argv[4])
qo_flag     = int(sys.argv[5])
n_processor = int(sys.argv[6])
cb_atoms_all=Zheng_func.get_cb_PDB_trj(pdb_file)

n_snap = len(cb_atoms_all)
print "Number of snapshots in the PDB file is: ", n_snap 

#constants
cutoff = 9.5
out = open(output_file,'w')

if n_processor > 1:
	import multiprocessing
	pool = multiprocessing.Pool(n_processor)

for i in range(n_snap):
	cb_atoms_1 = cb_atoms_all[i]
	
	lists_1_2 = []
	
	for j in range(i+1,n_snap):
		cb_atoms_2 = cb_atoms_all[j]
	
		lists_1_2.append([cb_atoms_1, cb_atoms_2])
	
		if n_processor == 1:					
			Q_ij,perm_list = Zheng_func.computeQ_permutation(cb_atoms_1, cb_atoms_2, n_chains, perm_flag, qo_flag, cutoff)
			out.write(str(round(1-Q_ij,5))+' ')	
	
	if n_processor > 1:
	##Use multiprocessing
		D_ij_list = pool.map(computeD_multi_processor, lists_1_2)		
		D_rounded = [ round(elem,5) for elem in D_ij_list ]	
		out.write(' '.join(map(str,D_rounded)))
		out.write(' ')
	
	sys.stdout.write(str(i)+' ')
	sys.stdout.flush()
out.write('\n')
out.close()


