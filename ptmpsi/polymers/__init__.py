import numpy as np
from ptmpsi.constants import covrad, ztosym, symtoz

def distance(a,b):
    return np.linalg.norm([float(x)-float(y) for x,y in zip(a,b)])

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
            dist = np.linalg.norm(ri-rj)
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

def rename(coords,bonds,natoms):
    for iatom in range(natoms):
        if coords[iatom,0] == 'C':
            obond = False
            nbond = False
            for jatom in range(natoms):
                if bonds[iatom,jatom]:
                    if coords[jatom,0] == 'O': obond = True
                    if coords[jatom,0] == 'N': nbond = True
            if nbond and not obond: coords[iatom,0] = 'CE'
            if nbond and obond: coords[iatom,0] = 'CO'
    for iatom in range(natoms):
        if coords[iatom,0] == 'C':
            cebond = False
            for jatom in range(natoms):
                if bonds[iatom,jatom]:
                    if coords[jatom,0] == 'CE': 
                        cebond = True
                        break
            if cebond: coords[iatom,0] = 'CD'
    for iatom in range(natoms):
        if coords[iatom,0] == 'C':
            cdbond = False
            for jatom in range(natoms):
                if bonds[iatom,jatom]:
                    if coords[jatom,0] == 'CD': 
                        cdbond = True
                        break
            if cdbond: coords[iatom,0] = 'CG'
    for iatom in range(natoms):
        if coords[iatom,0] == 'C':
            cgbond = False
            for jatom in range(natoms):
                if bonds[iatom,jatom]:
                    if coords[jatom,0] == 'CG': 
                        cgbond = True
                        break
            if cgbond: coords[iatom,0] = 'CB'
    for iatom in range(natoms):
        if coords[iatom,0] == 'C':
            cbbond = False
            for jatom in range(natoms):
                if bonds[iatom,jatom]:
                    if coords[jatom,0] == 'CB': 
                        cbbond = True
                        break
            if cbbond: coords[iatom,0] = 'CA'
    for iatom in range(natoms):
        if coords[iatom,0] == 'C':
            cbond = False
            for jatom in range(natoms):
                if bonds[iatom,jatom]:
                    if coords[jatom,0] == 'CO': 
                        cbond = True
                        break
            if cbond: coords[iatom,0] = 'CA'
    for iatom in range(natoms):
        if coords[iatom,0] == 'C':
            cbond = False
            for jatom in range(natoms):
                if bonds[iatom,jatom]:
                    if coords[jatom,0] == 'CA': 
                        cbond = True
                        break
            if cbond: coords[iatom,0] = 'CB'
    for iatom in range(natoms):
        if coords[iatom,0] == 'C':
            cbond = False
            for jatom in range(natoms):
                if bonds[iatom,jatom]:
                    if coords[jatom,0] == 'CB': 
                        cbond = True
                        break
            if cbond: coords[iatom,0] = 'CG'
    for iatom in range(natoms):
        if coords[iatom,0] == 'C':
            cbond = False
            for jatom in range(natoms):
                if bonds[iatom,jatom]:
                    if coords[jatom,0] == 'CG': 
                        cbond = True
                        break
            if cbond: coords[iatom,0] = 'CD'
    for iatom in range(natoms):
        if coords[iatom,0] == 'C':
            cbond = False
            for jatom in range(natoms):
                if bonds[iatom,jatom]:
                    if coords[jatom,0] == 'CD': 
                        cbond = True
                        break
            if cbond: coords[iatom,0] = 'CE'
    for iatom in range(natoms):
        if coords[iatom,0] == 'CO':
            coords[iatom,0] = 'C'
            
alphabet = {"CA":6,"CB":5, "CG":4, "CD":3, "CE":2, "N":1, "C":0}

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
        
dict2 = {"N": "CE", "CE":"CD", "CD":"CG", "CG":"CB", "CB":"CA", "CA":"C", "C":"O", "O":"N"}

def sortcoords(coords,nterminal,natoms,bonds):
    newcoords = []
    newcoords.append(coords[nterminal])
    done = np.empty(natoms,dtype=bool)
    done[:] = False
    done[nterminal] = True
    lastpos = nterminal
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
                    if bonded_frag[lastpos,iatom]:
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
        if newcoords[-1][0] == "CE":
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
def gen_xyz(rangea=1,rangeb=1,rangec=1):
    xyzcoords = []
    labels = []
    for a in range(rangea):
        for b in range(rangeb):
            for c in range(rangec):
                for coord in nylon6_coords:
                    vector1 = np.array([float(coord[1]),float(coord[2]),float(coord[3])]) 
                    vector2 = np.array([-vector1[0],-vector1[1],vector1[2]+0.5])
                    vector2 += np.array([a,b,c])
                    vector1 += np.array([a,b,c])
                    cart1 = Minv.T @ vector1
                    cart2 = Minv.T @ vector2
                    labels.append(coord[0])
                    labels.append(coord[0])
                    xyzcoords.append([cart1[0],cart1[1],cart1[2]])
                    xyzcoords.append([cart2[0],cart2[1],cart2[2]])
    temp = np.array(xyzcoords)
    templabels = np.array(labels)
    maxy = np.max(temp,axis=0)[1]
    miny = np.min(temp,axis=0)[1]
  
            
    for i,coord in enumerate(temp):
        if np.isclose(coord[1],miny,atol=0.3):
            xyzcoords.append(coord + rangeb*_lattice_vector_b)
            labels.append(labels[i])
        if np.isclose(coord[1],maxy,atol=0.3):
            xyzcoords.append(coord - rangeb*_lattice_vector_b)
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
                xyzcoords.append(_coord - rangea*_lattice_vector_a)
            if np.isclose(_coord[0],minx,atol=2.5):
                labels.append(label[i])
                xyzcoords.append(_coord + rangea*_lattice_vector_a)
                
    return xyzcoords,labels

# Data from European Polymer Journal 1979, 15, 765.
_lattice_constant_a = 9.71
_lattice_constant_b = 8.19
_lattice_constant_c = 17.24
_angle_alpha = np.radians(90.0)
_angle_beta = np.radians(90.0)
_angle_gamma = np.radians(112.5)

nylon6_coords = np.array([
   ["C", -0.052900,-0.003500,0.000000],
   ["C",  0.047000, 0.003200,0.069700],
   ["C", -0.036600,-0.002500,0.145600],
   ["C",  0.059200, 0.004000,0.214900],
   ["C",  0.069200, 0.004600,0.354400],
   ["C", -0.030700,-0.002000,0.424200],
   ["N ", -0.008200,-0.000600,0.281300],
   ["O ",  0.189300, 0.012700,0.208700],
   ["C",  0.447100,-0.003500,0.490000],
   ["C",  0.547000, 0.003200,0.420200],
   ["C",  0.463400,-0.002500,0.344400],
   ["C",  0.559200, 0.004000,0.275100],
   ["C",  0.569200, 0.004600,0.135500],
   ["C",  0.469300,-0.002000,0.065800],
   ["N ",  0.491800,-0.000600,0.208700],
   ["O ",  0.689300, 0.012700,0.281300],
   ["C", -0.052900, 0.496500,0.214400],
   ["C",  0.047000, 0.503100,0.284100],
   ["C", -0.036600, 0.497600,0.360000],
   ["C",  0.059200, 0.504000,0.429300],
   ["C",  0.069200, 0.504600,0.568800],
   ["C", -0.030700, 0.497900,0.638500],
   ["N ", -0.008200, 0.499400,0.495600],
   ["O ",  0.189300, 0.512700,0.423100],
   ["C",  0.447100, 0.496500,0.704300],
   ["C",  0.547000, 0.503100,0.634600],
   ["C",  0.463400, 0.497600,0.558700],
   ["C",  0.559200, 0.504000,0.489400],
   ["C",  0.569200, 0.504600,0.349900],
   ["C",  0.469300, 0.497900,0.280200],
   ["N ",  0.491800, 0.499400,0.423100],
   ["O ",  0.689300, 0.512700,0.4956]
])


# Define Lattice vectors
V = _lattice_constant_a * _lattice_constant_b * _lattice_constant_c * np.sqrt(1+2*np.cos(_angle_alpha)*np.cos(_angle_beta)*np.cos(_angle_gamma)-np.cos(_angle_alpha)**2-np.cos(_angle_beta)**2-np.cos(_angle_gamma)**2 )
_lattice_vector_a = _lattice_constant_a * np.array([1.0,0,0])
_lattice_vector_b = _lattice_constant_b * np.array([np.cos(_angle_gamma),np.sin(_angle_gamma),0.0])
_lattice_vector_c = _lattice_constant_c * np.array([np.cos(_angle_beta),(np.cos(_angle_alpha)-np.cos(_angle_beta)*np.cos(_angle_gamma))/np.sin(_angle_gamma),V/(_lattice_constant_a*_lattice_constant_b*_lattice_constant_c*np.sin(_angle_gamma))])
Minv = np.array([_lattice_vector_a,_lattice_vector_b,_lattice_vector_c])

# Function to generate XYZ coordinates
def gen_xyz(rangea=1,rangeb=1,rangec=1):
    xyzcoords = []
    labels = []
    for a in range(rangea):
        for b in range(rangeb):
            for c in range(rangec):
                for coord in nylon6_coords:
                    vector1 = np.array([float(coord[1]),float(coord[2]),float(coord[3])]) 
                    vector2 = np.array([-vector1[0],-vector1[1],vector1[2]+0.5])
                    vector2 += np.array([a,b,c])
                    vector1 += np.array([a,b,c])
                    cart1 = Minv.T @ vector1
                    cart2 = Minv.T @ vector2
                    labels.append(coord[0])
                    labels.append(coord[0])
                    xyzcoords.append([cart1[0],cart1[1],cart1[2]])
                    xyzcoords.append([cart2[0],cart2[1],cart2[2]])
    temp = np.array(xyzcoords)
    templabels = np.array(labels)
    maxy = np.max(temp,axis=0)[1]
    miny = np.min(temp,axis=0)[1]
  
            
    for i,coord in enumerate(temp):
        if np.isclose(coord[1],miny,atol=0.3):
            xyzcoords.append(coord + rangeb*_lattice_vector_b)
            labels.append(labels[i])
        if np.isclose(coord[1],maxy,atol=0.3):
            xyzcoords.append(coord - rangeb*_lattice_vector_b)
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
                xyzcoords.append(_coord - rangea*_lattice_vector_a)
            if np.isclose(_coord[0],minx,atol=2.5):
                labels.append(label[i])
                xyzcoords.append(_coord + rangea*_lattice_vector_a)
                
    return xyzcoords,labels

def gen_pdb(rangea=1, rangeb=1, rangec=1):
    xyzcoords, labels = gen_xyz(rangea,rangeb,rangec)

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
        coords_frag = sortcoords(coords_frag,nterminal,natoms_frag,bonded_frag)
        idatom = addchain(dict3[idfrag%62],idatom,coords_frag,"nylon6_{rangea}{rangeb}{rangec}.pdb")
        idatom = addchainmdl(idfrag+1,idatom,coords_frag,f"nylon6_models_{rangea}{rangeb}{rangec}.pdb")
    return
