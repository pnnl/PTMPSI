# ----------------------------------------------------------------------
# Copyright (2013) Weihua Zheng@Wolynes group, Rice University
# Last Update: 03/16/2013
# ----------------------------------------------------------------------
from Bio.PDB.PDBParser import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
from math import *
import numpy ## array_2d =  zeros([10 10], float)
import os, sys
from itertools import permutations
from itertools import combinations
from scipy.spatial.distance import pdist, squareform
#from MDAnalysis.core.util import norm

atom_desc = {'1' : 'C-Alpha', '2' : 'N', '3' : 'O', '4' : 'C-Beta', '5' : 'H-Beta', '6' : 'C-Prime'}
PDB_type = {'1' : 'CA', '2' : 'N', '3' : 'O', '4' : 'CB', '5' : 'HB', '6' : 'C' }

def vector(p1, p2):
    return [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]

def vabs(a):
    #return sqrt(pow(a[0],2)+pow(a[1],2)+pow(a[2],2)) #much slower
    return sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] )

def sigma_sq(sep):
	return pow((1+sep),0.15)*pow((1+sep), 0.15)

def calc_dist(p1, p2):
	v=vector(p1,p2)
	return vabs(v)

def save2d(filename, array_2d, format):
	numpy.savetxt(filename, array_2d, fmt=format);

def get_file_len(filename):
    exeline="wc -l "+filename ; 
    stdout=os.popen(exeline).readline(); 
    stdout=stdout.split()
    return int(stdout[0])

def getline_range(filename, line1, line2):
    assert(line1 <= line2)
    nline=line2-line1+1
    exeline="head -n"+str(line2)+" "+filename+" | tail -n"+str(nline) ; 
    #print exeline
    stdout=os.popen(exeline).readlines()
    return stdout
def get_natoms(dump_file):
	line=getline_range(dump_file, 4,4); line=line[0].split()
	natoms=int(line[0])
	return natoms
def get_nline_snapshot(dump_file):
	natoms=get_natoms(dump_file)
	nline=natoms + 9
	return nline

def get_dump_i(dump_file, i):
	nline=get_nline_snapshot(dump_file)
	line_start = 1 + nline*i ; line_end = nline*(i+1)
	dump_part=getline_range(dump_file, line_start, line_end)
	return dump_part

def pdbname(argv_i):
	name=argv_i
	if name[-4:].lower()==".pdb":
		pdb=name
	else:
		pdb=name+".pdb"
	return pdb

def compute_N_contacts(cb_atoms, cutoff, min_sep):
        count = 0
        N = len(cb_atoms)
        for ia in range(0, N):
                for ja in range(ia+min_sep, N):
                        r = vabs(vector(cb_atoms[ia], cb_atoms[ja]))
			if r <= cutoff: count += 1
        return count

def compute_inter_contacts(cb_atoms, chainIDs, cutoff, min_sep):
        count = 0
        N = len(cb_atoms)
        for ia in range(0, N):
                for ja in range(ia+min_sep, N):
			if chainIDs[ia] != chainIDs[ja]:
                        	r = vabs(vector(cb_atoms[ia], cb_atoms[ja]))
				if r <= cutoff : count += 1
        return count

def inter_contact_check(cb_atoms, ires_list, chainIDs, cutoff, min_sep):
        N   = len(cb_atoms)
	Nres= len(ires_list)
	ires_list_flag = [0]*Nres
        for i in range(Nres):
		ia = ires_list[i]
		ja = 0
		while ires_list_flag[i] == 0 and ja < N:
			if chainIDs[ia] != chainIDs[ja]:
                        	r = vabs(vector(cb_atoms[ia], cb_atoms[ja]))
				if r <= cutoff : ires_list_flag[i] = 1
			ja = ja + 1
	return ires_list_flag

def compute_contact_order(cb_atoms, cutoff, min_sep):
        count = 0; contact_order = 0
        N = len(cb_atoms)
        for ia in range(0, N):
                for ja in range(ia+min_sep, N):
                        r = vabs(vector(cb_atoms[ia], cb_atoms[ja]))
			if r <= cutoff: 
				count += 1
				contact_order += abs(ia-ja)
	count=max(count,1) #to avoid it being zero
	return contact_order / float(count)

def get_n_chains(pdbfile):
	p=PDBParser(PERMISSIVE=1)
	s = p.get_structure(pdbfile,pdbfile)
	chains = s[0].get_list()
	return len(chains)

def get_chainIDs_for_residues(pdbfile) :
	p=PDBParser(PERMISSIVE=1)
	s = p.get_structure(pdbfile,pdbfile)
	chains = s[0].get_list()
	chainIDs = []
	i_chain=65
	for chain in chains:
		ID = chain.get_id()
		for res in chain:
			chainIDs.append(ID)
	return chainIDs

def get_ca(pdbfile):
	p=PDBParser(PERMISSIVE=1)
	ca_atoms = []
	s = p.get_structure(pdbfile,pdbfile)
	chains = s[0].get_list()
	for chain in chains:
        	for res in chain:
               		is_regular_res = res.has_id('CA') and res.has_id('O')
	                res_id = res.get_id()[0]
        	        if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS' ) and is_regular_res:
                	        resname = res.get_resname(); 
                                ca_atoms.append(res['CA'].get_coord())
                	else :
                        	print "Pdb file contains irregular residue names or missing CA / O atoms! Fix it and run again! Exit with error."
				print "res_id :", res_id
				sys.exit()
	return ca_atoms

def get_cb_PDB_trj(pdb_trj):
    fh=open(pdb_trj, 'r')
    tmp_filename = 'tmp135791113.pdb'
    tmp_fh=open(tmp_filename, 'w')
    cb_atoms_all=[]
    #select residues to create cb_atoms lists
    for line in fh.readlines():
        head = line[0:3]    
        tmp_fh.write(line)
        
        if head == 'END':
            tmp_fh.close()
            cb_atoms = get_cbeta(tmp_filename)        
            cb_atoms_all.append(cb_atoms)
            tmp_fh=open(tmp_filename, 'w')
    tmp_fh.close()
    return cb_atoms_all
    
def get_cb(pdbfile):  ##same as get_cbeta
	p=PDBParser(PERMISSIVE=1)
	cb_atoms = []
	s = p.get_structure(pdbfile,pdbfile)
	chains = s[0].get_list()
	for chain in chains:
        	for res in chain:
               		is_regular_res = res.has_id('CA') and res.has_id('O')
	                res_id = res.get_id()[0]
        	        if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS' ) and is_regular_res:
                	        resname = res.get_resname(); 
                                if resname == 'GLY': cb_atoms.append(res['CA'].get_coord())
                                else :               cb_atoms.append(res['CB'].get_coord())
                	else :
                        	print "Pdb file contains irregular residue names or missing CA / O atoms! Fix it and run again! Exit with error."
				print "res_id :", res_id
				sys.exit()
	return cb_atoms

def get_cbeta(pdbfile):
	p=PDBParser(PERMISSIVE=1)
	cb_atoms = []
	s = p.get_structure(pdbfile,pdbfile)
	chains = s[0].get_list()
	for chain in chains:
        	for res in chain:
               		is_regular_res = res.has_id('CA') and res.has_id('O')
	                res_id = res.get_id()[0]
        	        if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS' ) and is_regular_res:
                	        resname = res.get_resname(); 
                                if resname == 'GLY': cb_atoms.append(res['CA'].get_coord())
                                else :               cb_atoms.append(res['CB'].get_coord())
                	else :
                        	print "Pdb file contains irregular residue names or missing CA / O atoms! Fix it and run again! Exit with error."
				print "res_id :", res_id
				sys.exit()
	return cb_atoms

def get_dump_n_snapshot(dump_file) :
	file_len=get_file_len(dump_file)
	line=getline_range(dump_file, 4,4); line=line[0].split()
	n_atoms = int(line[0]); nline=n_atoms + 9
	n_snapshot=file_len / nline
	return n_snapshot, nline

def get_dump_boxsize_x(dump_file) :
	#file_len=get_file_len(dump_file)
	line=getline_range(dump_file, 6,6); line=line[0].split()
	x_1 = float(line[0]); 
	x_2 = float(line[1]); 
	dx = x_2 - x_1
	return dx

def get_atoms_dump(dump_lines):
	ca_atoms_dump = []; box = []; A = []
	cb_atoms_dump = [] 
	o_atoms_dump = [] 
	for l in dump_lines:
		l = l.strip()
		if l[:5]=="ITEM:": item = l[6:]
		elif item[:10] == "BOX BOUNDS":
		    box.append(l)
		    l = l.split()
		    A.append([float(l[0]), float(l[1])])
		elif item[:5] == "ATOMS":
			l = l.split()  ; i_atom = l[0]
			x = float(l[2]); y = float(l[3]); z = float(l[4])
			x = (A[0][1] - A[0][0])*x + A[0][0]
			y = (A[1][1] - A[1][0])*y + A[1][0]
			z = (A[2][1] - A[2][0])*z + A[2][0]
			desc = atom_desc[l[1]]
			if desc=='C-Alpha':
				atom_CA = [x,y,z]
				ca_atoms_dump.append(atom_CA)
			if desc=='C-Beta':
				atom = [x,y,z]
				cb_atoms_dump.append(atom) 
			if desc=='H-Beta' :
				cb_atoms_dump.append(atom_CA)
			if desc=='O' :
				atom = [x,y,z]
				o_atoms_dump.append(atom)
	return ca_atoms_dump, cb_atoms_dump, o_atoms_dump

def get_ca_dump(dump_lines):
	ca_atoms_dump = []; box = []; A = []
	for l in dump_lines:
		l = l.strip()
		if l[:5]=="ITEM:": item = l[6:]
		elif item[:10] == "BOX BOUNDS":
		    box.append(l)
		    l = l.split()
		    A.append([float(l[0]), float(l[1])])
		elif item[:5] == "ATOMS":
			l = l.split()  ; i_atom = l[0]
			x = float(l[2]); y = float(l[3]); z = float(l[4])
			x = (A[0][1] - A[0][0])*x + A[0][0]
			y = (A[1][1] - A[1][0])*y + A[1][0]
			z = (A[2][1] - A[2][0])*z + A[2][0]
			desc = atom_desc[l[1]]
			if desc=='C-Alpha':
				atom_CA = [x,y,z]
				ca_atoms_dump.append(atom_CA)
	return ca_atoms_dump

def get_cbeta_dump(dump_lines):
	cb_atoms_dump = []; box = []; A = []
	for l in dump_lines:
		l = l.strip()
		if l[:5]=="ITEM:": item = l[6:]
		elif item[:10] == "BOX BOUNDS":
		    box.append(l)
		    l = l.split()
		    A.append([float(l[0]), float(l[1])])
		elif item[:5] == "ATOMS":
			l = l.split()  ; i_atom = l[0]
			x = float(l[2]); y = float(l[3]); z = float(l[4])
			x = (A[0][1] - A[0][0])*x + A[0][0]
			y = (A[1][1] - A[1][0])*y + A[1][0]
			z = (A[2][1] - A[2][0])*z + A[2][0]
			desc = atom_desc[l[1]]
			if desc=='C-Alpha':
				atom_CA = [x,y,z]
				#ca_atoms.append(atom)
			if desc=='C-Beta':
				atom = [x,y,z]
				cb_atoms_dump.append(atom) 
			if desc=='H-Beta' :
				cb_atoms_dump.append(atom_CA)
	return cb_atoms_dump
def fused_count_contact(cb_atoms_pdb, cb_atoms, Nmer, len_linker, cutoff):
	Ntotal     = len(cb_atoms); 
	N_mono = len(cb_atoms_pdb); 
	N_A = (Ntotal - (Nmer-1)*len_linker)/Nmer
	if N_mono != N_A :
		print "N_mono=", N_mono, "N_A=", N_A
		print "Only apply this code to domains with the same length! Exit with error!"
		sys.exit()
	#print Ntotal, N_mono, Nmer, len_linker

	l_native=[0]*Nmer*2  #native, intra-non-native 
	l_nonnative=numpy.zeros([Nmer*(Nmer-1)/2, 4], int) #2d array, 4d: domain-swap, self, self-shifted, other-inter
	min_sep = 4
	count=0
	for i in range(Nmer):
		#start_shift=max(0, res_start-1); end_shift=min(N_mono, res_end)
		i_start = i*(N_mono+len_linker); i_end = i_start + N_mono
		#i_start = i*(N_mono+len_linker)+start_shift ; i_end = i_start + end_shift
		for j in range(i,Nmer):
			j_start = j*(N_mono+len_linker); j_end = j_start + N_mono 
			for ia in range(i_start, i_end):
				for ja in range(j_start,j_end):
					if i==j :
						if ja-ia < min_sep : continue
						else :
							r1=  vabs(vector(cb_atoms[ia],     cb_atoms[ja]))
							if r1 <= cutoff:
								ia_mono = ia - i_start ; ja_mono = ja - j_start 
								rn = vabs(vector(cb_atoms_pdb[ia_mono], cb_atoms_pdb[ja_mono]))
								if rn <= cutoff : l_native[i]+=1   #native
								if rn >  cutoff : l_native[i+Nmer]+=1   #intra-nonnative
					else :
						r1=  vabs(vector(cb_atoms[ia],     cb_atoms[ja]))
						if r1 <= cutoff:
							ia_mono = ia - i_start ; ja_mono = ja - j_start 
							if ia_mono==ja_mono: 
								l_nonnative[count][1]+=1 #self-recognized
							elif fabs(ia_mono - ja_mono) < min_sep:
								l_nonnative[count][2]+=1 #self shifted inter-contacts
							elif fabs(ia_mono - ja_mono) >= min_sep:
								rn = vabs(vector(cb_atoms_pdb[ia_mono], cb_atoms_pdb[ja_mono]))
								if rn > cutoff : l_nonnative[count][3]+=1 #other inter-contacts
								else : l_nonnative[count][0]+=1 #domain_swap
			if i != j : count += 1
	#return l_mono, l_inter
	return l_native, l_nonnative
###For homogenous Oligomer systems, counting the number of intra/inter-monomer beta HBonds
def count_beta_contact(ca_atoms_pdb, ca_atoms, Nmer, len_linker, cutoff):
	Ntotal     = len(ca_atoms); 
	N_mono = len(ca_atoms_pdb); 
	N_A = (Ntotal - (Nmer-1)*len_linker)/Nmer
	if N_mono != N_A :
		print "N_mono=", N_mono, "N_A=", N_A
		print "Homogenous peptides are required!"
		print "Only apply this code to domains with the same length! Exit with error!"
		sys.exit()
	#print Ntotal, N_mono, Nmer, len_linker

	l_intra=[0]*2  #intra-antipara, intra-para 
	l_inter=[0]*2  #inter-antipara, inter-para 
	min_sep = 5
	count=0
	for i in range(Nmer):
		#start_shift=max(0, res_start-1); end_shift=min(N_mono, res_end)
		i_start = i*(N_mono+len_linker); i_end = i_start + N_mono
		#i_start = i*(N_mono+len_linker)+start_shift ; i_end = i_start + end_shift
		for j in range(i,Nmer):
			j_start = j*(N_mono+len_linker); j_end = j_start + N_mono 
			for ia in range(i_start, i_end):
				for ja in range(j_start,j_end):
					if i==j : #intra
						if ja-ia < min_sep : continue
						else :
							ia_mono = ia - i_start ; ja_mono = ja - j_start 
							if ia_mono == 0 or ja_mono == N_mono -1 : continue
							r1=  vabs(vector(ca_atoms[ia],     ca_atoms[ja]))
							r2=  vabs(vector(ca_atoms[ia-1],   ca_atoms[ja+1]))
							r3=  vabs(vector(ca_atoms[ia+1],   ca_atoms[ja+1]))
							if r1 <= cutoff and r2 <= cutoff : #anti-para
								l_intra[0]+=1  
							if r1 <= cutoff and r3 <= cutoff : #para
								l_intra[1]+=1  
					else :   #inter
						ia_mono = ia - i_start ; ja_mono = ja - j_start 
						if ia_mono == 0 or ja_mono == N_mono -1 : continue
						r1=  vabs(vector(ca_atoms[ia],     ca_atoms[ja]))
						r2=  vabs(vector(ca_atoms[ia-1],   ca_atoms[ja+1]))
						r3=  vabs(vector(ca_atoms[ia+1],   ca_atoms[ja+1]))
						if r1 <= cutoff and r2 <= cutoff : #anti-para
							l_inter[0]+=1  
						if r1 <= cutoff and r3 <= cutoff : #para
							l_inter[1]+=1
	return l_intra, l_inter

def fused_count_contact_self(cb_atoms_pdb, cb_atoms, Nmer, len_linker, cutoff, res_start, res_end):
	Ntotal     = len(cb_atoms); 
	N_mono = len(cb_atoms_pdb); 
	N_A = (Ntotal - (Nmer-1)*len_linker)/Nmer
	if N_mono != N_A :
		print "N_mono=", N_mono, "N_A=", N_A
		print "Only apply this code to domains with the same length! Exit with error!"
		sys.exit()
	#print Ntotal, N_mono, Nmer, len_linker

	l_native=[0]*Nmer*2  #native, intra-non-native 
	l_nonnative=numpy.zeros([Nmer*(Nmer-1)/2, 4], int) #2d array, 4d: domain-swap, self, self-shifted, other-inter
	min_sep = 4
	count=0
	for i in range(Nmer):
		start_shift=max(0, res_start-1); end_shift=min(N_mono, res_end)
		#i_start = i*(N_mono+len_linker); i_end = i_start + N_mono
		i_start = i*(N_mono+len_linker)+start_shift ; i_end = i*(N_mono+len_linker) + end_shift
		for j in range(i,Nmer):
			j_start = j*(N_mono+len_linker)+start_shift; j_end = j*(N_mono+len_linker) + end_shift
			for ia in range(i_start, i_end):
				for ja in range(j_start,j_end):
					if i==j :
						if ja-ia < min_sep : continue
						else :
							r1=  vabs(vector(cb_atoms[ia],     cb_atoms[ja]))
							if r1 <= cutoff:
								ia_mono = ia - i_start ; ja_mono = ja - j_start 
								rn = vabs(vector(cb_atoms_pdb[ia_mono], cb_atoms_pdb[ja_mono]))
								if rn <= cutoff : l_native[i]+=1   #native
								if rn >  cutoff : l_native[i+Nmer]+=1   #intra-nonnative
					else :
						#print start_shift, end_shift, i_start, i_end, j_start, j_end, ia, ja; 
						r1=  vabs(vector(cb_atoms[ia],     cb_atoms[ja])); 
						if r1 <= cutoff:
							ia_mono = ia - i_start ; ja_mono = ja - j_start 
							if ia_mono==ja_mono: 
								l_nonnative[count][1]+=1 #self-recognized
							elif fabs(ia_mono - ja_mono) < min_sep:
								l_nonnative[count][2]+=1 #self shifted inter-contacts
							elif fabs(ia_mono - ja_mono) >= min_sep:
								rn = vabs(vector(cb_atoms_pdb[ia_mono], cb_atoms_pdb[ja_mono]))
								if rn > cutoff : l_nonnative[count][3]+=1 #other inter-contacts
								else : l_nonnative[count][0]+=1 #domain_swap
			if i != j : count += 1
	#return l_mono, l_inter
	return l_native, l_nonnative

def contactmap(cb_atoms, cutoff):
	#c_native = 0.25 ; c_boundary=0.1 ; c_both   = 1.0 ; c_dump   = 0.6
	Ntotal     = len(cb_atoms); 
	#N_mono = len(cb_atoms_pdb); 
	#N_A = (Ntotal - (Nmer-1)*len_linker)/Nmer
	#if N_mono != N_A :
	#	print "N_mono=", N_mono, "N_A=", N_A
	#	print "Only apply this code to domains with the same length! Exit with error!"
	#	sys.exit()
	contact_matrix = numpy.zeros([Ntotal, Ntotal], float)
	min_sep = 3
	for i in range(Ntotal):
		for j in range(i+min_sep,Ntotal):
			r=  vabs(vector(cb_atoms[i], cb_atoms[j]))
			if r <= cutoff : contact_matrix[i, j] = 1
	return contact_matrix

def fused_contactmap(cb_atoms_pdb, cb_atoms, Nmer, len_linker, cutoff):
	c_native = 0.25 ; c_boundary=0.1 ; c_both   = 1.0 ; c_dump   = 0.6

	Ntotal     = len(cb_atoms); 
	N_mono = len(cb_atoms_pdb); 
	N_A = (Ntotal - (Nmer-1)*len_linker)/Nmer
	if N_mono != N_A :
		print "N_mono=", N_mono, "N_A=", N_A
		print "Only apply this code to domains with the same length! Exit with error!"
		sys.exit()
	contact_matrix = numpy.zeros([Ntotal, Ntotal], float)
	min_sep = 3
	#linkers are excluded when forming contact_map
	for i in range(Nmer):
		i_start = i*(N_mono+len_linker) ; i_end = i_start + N_mono 
		for j in range(i,Nmer):
			j_start = j*(N_mono+len_linker) ; j_end = j_start + N_mono 
			for ia in range(i_start, i_end):
				for ja in range(j_start,j_end):
					if i==j :
						if ja-ia < min_sep : continue
						else :
							ia_mono = ia - i_start ; ja_mono = ja - j_start 
							rn = vabs(vector(cb_atoms_pdb[ia_mono], cb_atoms_pdb[ja_mono]))
							if rn <= cutoff :
								contact_matrix[ja, ia] = c_both #record native contact map in the upper diagonal
							r1=  vabs(vector(cb_atoms[ia],     cb_atoms[ja]))
							if r1 <= cutoff and rn <= cutoff : contact_matrix[ia, ja] = c_both #l_native[i]+=1   #native
							if r1 <= cutoff and rn >  cutoff : contact_matrix[ia, ja] = c_dump 
							if r1 >  cutoff and rn <= cutoff : contact_matrix[ia, ja] = c_native 
					else :
						r1=  vabs(vector(cb_atoms[ia],     cb_atoms[ja]))
						if r1 <= cutoff: contact_matrix[ia, ja] = c_dump
						ia_mono = ia - i_start ; ja_mono = ja - j_start 
						rn = vabs(vector(cb_atoms_pdb[ia_mono], cb_atoms_pdb[ja_mono]))
						if rn <= cutoff and fabs(ja_mono-ia_mono) >= min_sep : 
							contact_matrix[ia, ja] = c_native #domain_swap
							if r1 <= cutoff : contact_matrix[ia, ja] = c_both
	for ia in range(Ntotal):
		for ja in range(Ntotal):
			ia_mono = ia%(N_A+len_linker)
			ja_mono = ja%(N_A+len_linker)
			if ia_mono > N_A -1 or ja_mono  > N_A -1 : contact_matrix[ia,ja]=c_boundary
			if ia_mono == ja_mono and contact_matrix[ia,ja]== 0 : contact_matrix[ia,ja]=c_boundary

	return contact_matrix

def fused_beta_contactmap(contact_matrix, beta_length):
        c_native = 0.25 ; c_boundary=0.1 ; c_both   = 1.0 ; c_dump   = 0.6
        Ntotal     = len(contact_matrix);
        betacontact_matrix = numpy.zeros([Ntotal, Ntotal], float)
        min_sep = 3
        for ia in range(Ntotal):
		betacontact_matrix[ia,ia] = c_boundary
                for ja in range(Ntotal):
                        if fabs(ja-ia) < min_sep : continue
                        if contact_matrix[ia, ja] == c_boundary : betacontact_matrix[ia,ja] = c_boundary
                        elif contact_matrix[ia,ja] >= c_dump :
                                #parallel
				count = 1; 
				ki =0; flag_contact = 1
				while (flag_contact) :
					ki-=1
					if ia+ki < 0 or ia+ki >= Ntotal or ja+ki < 0 or ja+ki >= Ntotal : flag_contact = 0
					elif contact_matrix[ia+ki,ja+ki] <  c_dump: flag_contact=0
					elif contact_matrix[ia+ki,ja+ki] >= c_dump: flag_contact = 1; count +=1;
				ki =0; flag_contact = 1
				while (flag_contact) :
					ki+=1
					if ia+ki < 0 or ia+ki >= Ntotal or ja+ki < 0 or ja+ki >= Ntotal : flag_contact = 0
					elif contact_matrix[ia+ki,ja+ki] <  c_dump: flag_contact=0
					elif contact_matrix[ia+ki,ja+ki] >= c_dump: flag_contact = 1; count +=1;
				if count >= beta_length : betacontact_matrix[ia,ja] = contact_matrix[ia,ja]

                                #anti-parallel
				count = 1; 
				ki =0; flag_contact = 1
				while (flag_contact) :
					ki-=1
					if ia+ki < 0 or ja-ki < 0 or ia+ki >= Ntotal or ja-ki >= Ntotal or fabs(ja-ki - (ia+ki)) < min_sep : flag_contact = 0
					elif contact_matrix[ia+ki,ja-ki] < c_dump: flag_contact=0
					elif contact_matrix[ia+ki,ja-ki] >=c_dump: flag_contact = 1; count +=1;
				ki =0; flag_contact = 1
				while (flag_contact) :
					ki+=1
					if ia+ki < 0 or ja-ki < 0 or ia+ki >= Ntotal or ja-ki >= Ntotal or fabs(ja-ki - (ia+ki)) < min_sep : flag_contact = 0
					elif contact_matrix[ia+ki,ja-ki] <  c_dump: flag_contact=0
					elif contact_matrix[ia+ki,ja-ki] >= c_dump: flag_contact = 1; count +=1;
				if count >= beta_length : betacontact_matrix[ia,ja] = contact_matrix[ia,ja]

        return betacontact_matrix

def compute_RMSD(ca_atoms_pdb, ca_atoms):
        if len(ca_atoms)!=len(ca_atoms_pdb):
                print "Error. Length mismatch!"
                print "dump file len: ", len(ca_atoms), ";  pdb file len: ", len(ca_atoms_pdb)
                sys.exit()

        l = len(ca_atoms)
        fixed_coord  = numpy.zeros((l, 3))
        moving_coord = numpy.zeros((l, 3))

        for i in range(0, l):
                fixed_coord[i]  = numpy.array ([ca_atoms_pdb[i][0], ca_atoms_pdb[i][1], ca_atoms_pdb[i][2]])
                moving_coord[i] = numpy.array ([ca_atoms[i][0], ca_atoms[i][1], ca_atoms[i][2]])
        sup = SVDSuperimposer()
        sup.set(fixed_coord, moving_coord)
        sup.run()
        rms = sup.get_rms()
        return rms
def computeQ_relative(cb_atoms_1, cb_atoms_2, qo_flag, cutoff): #optimized 5 star code
    N  = len(cb_atoms_1)
    N2 = len(cb_atoms_2) 
    if N != N2 :
        print "Error. Length mismatch!"
        print "file 1 Nres: ", N, ";  file 2 Nres: ", N2
        sys.exit()
                
    if qo_flag == 1 : min_sep = 4
    else:             min_sep = 3
    
    r1_matrix = squareform(pdist(cb_atoms_1)) #scipy pair-dist function
    r2_matrix = squareform(pdist(cb_atoms_2)) #fast matrix operation 
        
    Q = 0; norm = 0;
    for ia in xrange(N):
        for ja in xrange(ia+min_sep, N):
            r1 = r1_matrix[ia][ja]
            r2 = r2_matrix[ia][ja]
            
            if (r1 < cutoff or r2 < cutoff) or qo_flag == 0 :		
                dr = r2 - r1
                qij = exp(-dr*dr/(2*sigma_sq(ja-ia)))
            
                Q    += qij
                norm += 1
                
    Q = Q / norm
    return Q               
            
def compute_RG(ca_atoms):
        N = len(ca_atoms)
	if N == 0: print "empty ca atom list. Exit with error."; sys.exit()
	RG = 0
	for ia in range(N):
		for ja in range(ia+1, N):
			rv = vector(ca_atoms[ia], ca_atoms[ja])
			rsq= pow(rv[0],2) + pow(rv[1],2) + pow(rv[2],2)
			RG = RG + rsq
	RG = sqrt(RG/N/N)
	return RG

def computeQdiff(ca_atoms, ca_atoms_pdb1, ca_atoms_pdb2):
	if len(ca_atoms)!=len(ca_atoms_pdb1):
		print "Error. Length mismatch!"
		print "Pdb1: ", len(ca_atoms_pdb1), "trj: ", len(ca_atoms)
		sys.exit()
	if len(ca_atoms)!=len(ca_atoms_pdb2):
		print "Error. Length mismatch!"
		print "Pdb2: ", len(ca_atoms_pdb2), "trj: ", len(ca_atoms)
		sys.exit()
	splitq = False
	Q = {}; norm = {}; N = len(ca_atoms); min_sep = 3
	if qo_flag == 1 : min_sep = 4
	for ia in range(0, N):
		for ja in range(ia+min_sep, N):
			if not splitq:
				rn1 = vabs(vector(ca_atoms_pdb1[ia], ca_atoms_pdb1[ja]))
				rn2 = vabs(vector(ca_atoms_pdb2[ia], ca_atoms_pdb2[ja]))
				if qo_flag == 1 and ( rn1 >= cutoff and rn2 >= cutoff ) : continue

				r = vabs(vector(ca_atoms[ia], ca_atoms[ja]))
				dr1 = r - rn1 ; dr2 = r - rn2
				if splitq: index = pdb_chain_id[ia]
				else: index = 1
				if not Q.has_key(index):
					Q[index] = 0.0
					norm[index] = 0
				Q[index] = Q[index] + exp(-dr1*dr1/(2*sigma_sq(ja-ia))) + 1 - exp(-dr2*dr2/(2*sigma_sq(ja-ia)))
				norm[index] = norm[index] + 1
	for key in Q:
		Q[key] = Q[key]/norm[key]
	return Q

def computeQ_permutation(cb_atoms_1, cb_atoms_2, n_chains, qo_flag, permutation_flag, cutoff):
    N = len(cb_atoms_1)
    cutoff=9.5    
    if N != len(cb_atoms_2):
        print "Error, cb_atoms length mismatch."
        print "Length of cb_atoms_1 is ", N, "Length of cb_atoms_2 is ", len(cb_atoms_2)
        print "Check the number of residues in the pdb files."
        sys.exit()
    if permutation_flag == 0:
        Q_rel = computeQ_relative(cb_atoms_1, cb_atoms_2, qo_flag, cutoff)
        perm_chainlist = range(n_chains)
        
    #create all permuated cb_atoms_2_perm list for Q_rel calculation, max Q_rel is selected.
    else: 
        N_perm = factorial(n_chains)
        Q_perm = []
        n_res = N / n_chains #number of res per chain
        assert(n_res == int(N/n_chains))
        
        #cool stuff
        N_reshape_list = [range(N)[i:i+n_res] for i in range(0,N,n_res)]
                
        n_chain_perm = list(permutations(range(n_chains)))
        
        for i_chainlist in n_chain_perm:
            perm_res_list = []
            for j in range(n_chains):
                i_chain = i_chainlist[j]
                perm_res_list.extend(N_reshape_list[i_chain])

            cb_atoms_2_perm = []
            for ia in range(N):
                cb_atoms_2_perm.append(cb_atoms_2[perm_res_list[ia]])

            Q = computeQ_relative(cb_atoms_1, cb_atoms_2_perm, qo_flag, cutoff)
            Q_perm.append(Q)
        Q_rel = max(Q_perm)
        id = Q_perm.index(Q_rel)
        perm_chainlist = n_chain_perm[id]
    return Q_rel, perm_chainlist

#Mod
#Define a function that compute Qrel for all subsets of n-x chains from a n-mer
#and choose the maximum Qrel.
def computeQ_subset(cb_atoms_1, k_chains, cb_atoms_2, n_chains, qo_flag, permutation_flag, cutoff):
    #total number of residues in each chain
    N_k = len(cb_atoms_1)
    N_n = len(cb_atoms_2)    
    if N_k > N_n :
        print "Error, first PDB has more residues than the second one. Switch them!"
        print "Length of PDB1 is ", N_k, "Length of PDB2 is ", N_n    
        sys.exit()

    else:
        n_res = N_k / k_chains
        assert(n_res == N_n/n_chains)
        
        combination_chain_lists       = list(combinations(range(n_chains),k_chains))
        N_n_reshape_list  = [range(N_n)[i:i+n_res] for i in range(0,N_n,n_res)]
        
        Q_sub               = []
        perm_chain_list_sub = []
        for i_chain_list in combination_chain_lists:
            subset_res_list   = []
             #convert chain list to residue list
            for j in i_chain_list:
                 subset_res_list.extend(N_n_reshape_list[j])
                 
            cb_atoms_2_subset = []
            for ia in range(N_k):
                cb_atoms_2_subset.append(cb_atoms_2[subset_res_list[ia]])            
             
            Q, perm_chain_list = computeQ_permutation(cb_atoms_1, cb_atoms_2_subset, k_chains, qo_flag, permutation_flag, cutoff)
            Q_sub.append(Q)
            perm_chain_list_sub.append(perm_chain_list)
            
        Q_rel      = max(Q_sub)
        id         = Q_sub.index(Q_rel)
        perm_order = perm_chain_list_sub[id]
        i_chain    = combination_chain_lists[id]
        
        perm_sub_chains = [i_chain[i] for i in perm_order]
        
    return Q_rel, perm_sub_chains

##include only intra-monomer contacts, exclude linker region and inter-monomer contacts
def fused_computeQ(ca_atoms_pdb, ca_atoms, Nmer, len_linker, cutoff):
	Ntotal = len(ca_atoms);
	Nmono  = len(ca_atoms_pdb)
	N = Nmer*Nmono+(Nmer-1)*len_linker
	if Ntotal != N :
		print "Number of residues don't match! Exit with error!"
		print "dump file: ", Ntotal, "calculated from mono pdb: ", N
		sys.exit()
	
	min_sep = 4; a=0; qsum=0;
	for ia in range(Ntotal):
		for ja in range(Ntotal):
			ia_mono=ia%(Nmono+len_linker);  ia_chain_id = ia / (Nmono+len_linker);
			ja_mono=ja%(Nmono+len_linker) ; ja_chain_id = ja / (Nmono+len_linker);
			if ia_chain_id == ja_chain_id :
				if (ia_mono < Nmono and ja_mono < Nmono ) and ja_mono >= ia_mono + min_sep :
					rN=calc_dist(ca_atoms_pdb[ia_mono], ca_atoms_pdb[ja_mono])
					if rN < cutoff: #only native contacts
						r =calc_dist(ca_atoms[ia]    , ca_atoms[ja])
						dr = r - rN
						q_ij = exp(-dr*dr/(2*sigma_sq(ja-ia)))
						qsum += q_ij
						a += 1.0
	a=max(1,a) #to avoid it being zero
	qsum = qsum / a
	return  qsum
def fused_contactmap_all(cb_atoms_pdb, cb_atoms, Nmer, len_linker, cutoff):
	c_native = 0.25 ; c_boundary=0.1 ; c_both   = 1.0 ; c_dump   = 0.6

	Ntotal     = len(cb_atoms); 
	N_mono = len(cb_atoms_pdb); 
	N_A = (Ntotal - (Nmer-1)*len_linker)/Nmer
	if N_mono != N_A :
		print "N_mono=", N_mono, "N_A=", N_A
		print "Only apply this code to domains with the same length! Exit with error!"
		sys.exit()
	contact_matrix = numpy.zeros([Ntotal, Ntotal], float)
	min_sep = 3
	#linkers are excluded when forming contact_map
	for i in range(Nmer):
		i_start = i*(N_mono+len_linker) ; i_end = i_start + N_mono 
		for j in range(Nmer):
			j_start = j*(N_mono+len_linker) ; j_end = j_start + N_mono 
			for ia in range(i_start, i_end):
				for ja in range(j_start,j_end):
					if i==j :
						if fabs(ja-ia) < min_sep : continue
						else :
							ia_mono = ia - i_start ; ja_mono = ja - j_start 
							rn = vabs(vector(cb_atoms_pdb[ia_mono], cb_atoms_pdb[ja_mono]))
							#if rn <= cutoff :
							#	contact_matrix[ja, ia] = c_both #record native contact map in the upper diagonal
							r1=  vabs(vector(cb_atoms[ia],     cb_atoms[ja]))
							if r1 <= cutoff and rn <= cutoff : contact_matrix[ia, ja] = c_both #l_native[i]+=1   #native
							if r1 <= cutoff and rn >  cutoff : contact_matrix[ia, ja] = c_dump
							if r1 >  cutoff and rn <= cutoff : contact_matrix[ia, ja] = c_native 
					else :
						r1=  vabs(vector(cb_atoms[ia],     cb_atoms[ja]))
						if r1 <= cutoff: contact_matrix[ia, ja] = c_dump
						ia_mono = ia - i_start ; ja_mono = ja - j_start 
						rn = vabs(vector(cb_atoms_pdb[ia_mono], cb_atoms_pdb[ja_mono]))
						if rn <= cutoff and fabs(ja_mono-ia_mono) >= min_sep : 
							contact_matrix[ia, ja] = c_native #domain_swap
							if r1 <= cutoff : contact_matrix[ia, ja] = c_both
	for ia in range(Ntotal):
		for ja in range(Ntotal):
			ia_mono = ia%(N_A+len_linker)
			ja_mono = ja%(N_A+len_linker)
			if ia_mono > N_A -1 or ja_mono  > N_A -1 : contact_matrix[ia,ja]=c_boundary
			if ia_mono == ja_mono and contact_matrix[ia,ja]== 0 : contact_matrix[ia,ja]=c_boundary

	return contact_matrix
