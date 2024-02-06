import numpy as np
from ptmpsi.constants import covrad, ztosym, symtoz

class UnitCell:
    def __init__(self, a, b, c, alpha, beta, gamma, coords, ref=None):
        self.a = a
        self.b = b
        self.c = c
        self.alpha = np.radians(alpha)
        self.beta = np.radians(beta)
        self.gamma = np.radians(gamma)
        self.coords = coords
        self.ref = ref

# Define Lattice vectors
        V = self.a * self.b * self.c * np.sqrt(1+2*np.cos(self.alpha)*np.cos(self.beta)*np.cos(self.gamma)-np.cos(self.alpha)**2-np.cos(self.beta)**2-np.cos(self.gamma)**2 )
        self.lattice_vector_a = self.a * np.array([1.0,0,0])
        self.lattice_vector_b = self.b * np.array([np.cos(self.gamma),np.sin(self.gamma),0.0])
        self.lattice_vector_c = self.c * np.array([np.cos(self.beta),
                         (np.cos(self.alpha)-np.cos(self.beta)*np.cos(self.gamma))/np.sin(self.gamma),
                         V/(self.a*self.b*self.c*np.sin(self.gamma))])
        self.Minv = np.array([self.lattice_vector_a,self.lattice_vector_b,self.lattice_vector_c])
        return


def distance(a,b):
    return np.sqrt(
      (float(a[0])-float(b[0]))**2 +
      (float(a[1])-float(b[1]))**2 +
      (float(a[2])-float(b[2]))**2)

def bondmat(natoms,coords):
    bonded = np.empty((natoms,natoms),dtype=bool)
    bonded[:,:] = False
    _temp = [ [coord[0], float(coord[1]), float(coord[2]), float(coord[3])] for coord in coords ]
    for iatom in range(natoms):
        bonds = 0
        icov = covrad[_temp[iatom][0]]
        skipat3 = True
        ri = np.array(_temp[iatom][1:4])
        print("Checking all bonds for atom: {}, out of: {}".format(iatom+1,natoms))
        for jatom in range(iatom,natoms):
            jcov = covrad[_temp[jatom][0]]
            rj = np.array(_temp[jatom][1:4])
            dist = distance(ri,rj)
            if  dist < 1.3*(icov+jcov):
                bonds += 1
                bonded[jatom,iatom] = True
                if _temp[jatom][0] in ['o', 'n']: skipat3 = False
                if skipat3 and bonds == 3: break
                if bonds == 4: break
    return bonded | bonded.T

def getfrag(natoms,bonded):
    nfrag = 0
    frag  = [0]*natoms
    frag[0] = 1
    totalat = 1
    while totalat < natoms:
        nfrag += 1
        finished = False
        while not finished:
            finished = True
            for iatom in range(natoms):
                if frag[iatom] != nfrag: continue
                for jatom in range(natoms):
                    if bonded[iatom,jatom] and frag[jatom] == 0:
                        finished = False
                        frag[jatom] = nfrag + 0
                        totalat += 1
        print("Completed strand number: {}".format(nfrag))
        for iatom in range(natoms):
            if frag[iatom] == 0:
                frag[iatom] = nfrag + 1
                totalat += 1
                break
    return np.array(frag)


def addchain(label,idatom,coords,pdbfile):
    with open(pdbfile,'a') as pdb:
        seqid = 1
        atomid = 0
        for atom in coords:
            if atom[0] == "O":
                elem = "O"
            elif atom[0] == "N":
                elem = "N"
                seqid += 1
            else:
                elem = "C"
            atomid += 1
            idatom += 1
            if seqid == 1:
                residue = "ACE"
            elif atomid >= len(coords) - 1:
                residue = "NME"
            else:
                residue = "EAC"
            
            pdb.write("ATOM   {0:>4} {1:>4} {7} {2}{8:>4}    {3:8.3f}{4:8.3f}{5:8.3f}                      {6:>2}\n".format(idatom,atom[0],label,float(atom[1]),float(atom[2]),float(atom[3]),elem,residue,seqid))
        pdb.write("TER\n")
        return idatom
        

def sortcoords(coords,nterminal,natoms,bonds,dict2):
    newcoords = []
    newcoords.append(coords[nterminal])
    done = np.empty(natoms,dtype=bool)
    done[:] = False
    done[nterminal] = True
    lastpos = nterminal
    lastname = "CE" if "CE" in dict2 else "CG"
    while len(newcoords) < natoms:
        find = dict2[newcoords[-1][0]]
        for iatom in range(natoms):
            if done[iatom]: continue
            if coords[iatom,0] == find:
                if find == "N":
                    if np.any(bonds[lastpos] & bonds[iatom]):
                        newcoords.append(coords[iatom])
                        done[iatom] = True
                        lastpos = iatom
                        break
                else:
                    if bonds[lastpos,iatom]:
                        newcoords.append(coords[iatom])
                        done[iatom] = True
                        lastpos = iatom
                        break
                        
    while True:
        if newcoords[0][0] == "CA":
            newcoords[0][0] = "CH3"
            break
        else:
            newcoords.pop(0)
    
    while True:
        if newcoords[-1][0] == lastname:
            newcoords[-1][0] = "CH3"
            break
        else:
            newcoords.pop(-1)
    return newcoords

dict3 = {i:chr(i+65) for i in range(26)}
for i in range(26):
    ii = i+26
    dict3[ii] = chr(i+97)
for i in range(10):
    ii = i+52
    dict3[ii] = chr(i+48)
    
def addchainmdl(mdlid,idatom,coords,pdbfile):
    with open(pdbfile,'a') as pdb:
        seqid = 1
        atomid = 0
        pdb.write("MODEL {}\n".format(mdlid))
        for atom in coords:
            if atom[0] == "O":
                elem = "O"
            elif atom[0] == "N":
                elem = "N"
                seqid += 1
            else:
                elem = "C"
            atomid += 1
            idatom += 1
            if seqid == 1:
                residue = "ACE"
            elif atomid >= len(coords) - 1:
                residue = "NME"
            else:
                residue = "EAC"

            pdb.write("ATOM   {0:>4} {1:>4} {7} {2}{8:>4}    {3:8.3f}{4:8.3f}{5:8.3f}                      {6:>2}\n".format(idatom,atom[0],"A",float(atom[1]),float(atom[2]),float(atom[3]),elem,residue,seqid))
        pdb.write("TER\n")
        pdb.write("ENDMDL\n")
        return idatom



# Function to generate XYZ coordinates
def gen_xyz(cell, rangea=1, rangeb=1, rangec=1, symmetry=True):
    xyzcoords = []
    labels = []

    # Define Lattice vectors
    for a in range(rangea):
        for b in range(rangeb):
            for c in range(rangec):
                for coord in cell.coords:
                    vector1 = np.array([float(coord[1]),float(coord[2]),float(coord[3])]) 
                    vector1 += np.array([a,b,c])
                    cart1 = cell.Minv.T @ vector1
                    labels.append(coord[0])
                    xyzcoords.append([cart1[0],cart1[1],cart1[2]])
                    if symmetry:
                        vector1 -= np.array([a,b,c])
                        vector2 = np.array([-vector1[0],-vector1[1],vector1[2]+0.5])
                        vector2 += np.array([a,b,c])
                        cart2 = cell.Minv.T @ vector2
                        labels.append(coord[0])
                        xyzcoords.append([cart2[0],cart2[1],cart2[2]])

    if not symmetry: return xyzcoords, labels

    temp = np.array(xyzcoords)
    templabels = np.array(labels)
    maxy = np.max(temp,axis=0)[1]
    miny = np.min(temp,axis=0)[1]
  
            
    for i,coord in enumerate(temp):
        if np.isclose(coord[1],miny,atol=0.3):
            xyzcoords.append(coord + rangeb*cell.lattice_vector_b)
            labels.append(labels[i])
        if np.isclose(coord[1],maxy,atol=0.3):
            xyzcoords.append(coord - rangeb*cell.lattice_vector_b)
            labels.append(labels[i])

    temp = np.array(xyzcoords)
    templabels = np.array(labels)
    nlayers = 1
    curry = temp[0,1]
    layer = [0] * len(temp)
    layer[0] = 1
    totat = 1
    while True:
        for i,coord in enumerate(temp):
            if layer[i] == 0:
                if np.isclose(coord[1],curry,atol=0.3):
                    layer[i] = nlayers
                    totat += 1
        if totat == len(temp): break
        for i,coord in enumerate(temp):
            if layer[i] == 0:
                curry = temp[i,1]
                nlayers += 1
                layer[i] = nlayers
                totat += 1
                break      
            
    for ilayer in range(nlayers):
        coord = temp[np.array(layer) == ilayer+1]
        label = templabels[np.array(layer) == ilayer+1]
        maxx = np.max(coord,axis=0)[0]
        minx = np.min(coord,axis=0)[0]
        for i,_coord in enumerate(coord):
            if np.isclose(_coord[0],maxx,atol=2.5):
                labels.append(label[i])
                xyzcoords.append(_coord - rangea*cell.lattice_vector_a)
            if np.isclose(_coord[0],minx,atol=2.5):
                labels.append(label[i])
                xyzcoords.append(_coord + rangea*cell.lattice_vector_a)
                
    return xyzcoords,labels


def gen_pdb(xyzcoords, labels, alphabet, dict2, rename, head, tail):
    # Generate PDB file
    natoms = len(xyzcoords)
    coords = np.array([[lab.strip().capitalize(),x[0],x[1],x[2]] for lab,x in zip(labels,xyzcoords)])
    bonded = bondmat(natoms,coords)
    fragment = getfrag(natoms,bonded)
    nfrag = np.max(fragment)
    idatom = 0
    for idfrag in range(nfrag):
        print("Generating strand: {}, out of: {}".format(idfrag+1,nfrag))
        coords_frag = coords[ fragment == idfrag+1 ]
        natoms_frag = len(coords_frag)
        bonded_frag = bondmat(natoms_frag,coords_frag)
        rename(coords_frag,bonded_frag,natoms_frag)
        
        ends = []
        for iatom in range(natoms_frag):
            bonds = 0
            for jatom in range(natoms_frag):
                if bonded_frag[iatom,jatom]: bonds += 1
            if bonds == 2 and coords_frag[iatom,0] != "O": ends.append(iatom)

        if head == "nylon4":
            if len(ends) == 1:
                if coords_frag[-1,0] != "O": 
                    ends.append(natoms_frag-1)
                else:
                    ends.append(natoms_frag-2)
    
        anum = alphabet[coords_frag[ends[0],0]]
        bnum = alphabet[coords_frag[ends[1],0]]
        nterminal = None
        find = dict2[coords_frag[ends[0],0]]
        for iatom in range(natoms_frag):
            if bonded_frag[ends[0],iatom]:
                if coords_frag[iatom,0] == find:
                    nterminal = ends[0]
                    break
        if nterminal is None: nterminal = ends[1]
        coords_frag = sortcoords(coords_frag,nterminal,natoms_frag,bonded_frag,dict2)
        idatom = addchain(dict3[idfrag%62],idatom,coords_frag,f"{head}_{tail}.pdb")
        idatom = addchainmdl(idfrag+1,idatom,coords_frag,f"{head}_models_{tail}.pdb")
    return
