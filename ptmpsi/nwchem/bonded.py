from ptmpsi.nwchem.templates import hessnw
from ptmpsi.nwchem.reader import readoptim, readhess
from ptmpsi.constants import ang2bohr, nm2bohr, eh2kjmol, covrad, mass
import numpy as np
import copy

def get_hessian(files=["alpha","beta"],charge=0,mult=1,memory=2000,geometry=None):
    """
    Prepares hessian calculation input files
    """
    for filename in files:
        if geometry is None:
            _geometry = "load "+filename+".xyz"
        with open(filename+"_hess.nw","w") as fh:
            fh.write(hessnw.format(
                memory=memory,
                name=filename,
                mult=mult,
                charge=charge,
                geometry=_geometry))
    return

def bondedks(filename,elements,geometry):
    """
    Computes bond and angle force constants using the modified Seminario method

    The output units are compatible with GROMACS, i.e.
        kj/(mol nm^2) or kj/(mol rad^2)
    """
    hessian,natoms = readhess(filename+"_hess.log")

    # NWChem hessian is in Eh/bohr^2/Kamu
    hessian *= 0.001 * eh2kjmol * ang2bohr**2 

    if len(elements) != natoms:
        print(" Error: number of atoms mismatch")
        quit()

    # Remove mass weighting
    for ixyz in range(3*natoms):
        for jxyz in range(3*natoms):
            hessian[ixyz,jxyz] *= np.sqrt(mass[elements[ixyz//3]] * mass[elements[jxyz//3]])

    # Get bond list
    bond_list, angle_list = bond_angle_list(elements,geometry)

    # Get eigenvalues and eigenvectors
    eigenvalues = np.zeros((natoms,natoms,3))
    eigenvectors = np.zeros((natoms,natoms,3,3))
    bond_lengths = np.zeros((natoms,natoms))
    for iatom in range(natoms):
        for jatom in range(natoms):
            bond_lengths[iatom,jatom] = np.linalg.norm(geometry[iatom] - geometry[jatom])
            matrix = np.ascontiguousarray(hessian[ iatom*3:(iatom+1)*3, jatom*3:(jatom+1)*3 ])
            eigenvalues[iatom,jatom],eigenvectors[iatom,jatom] = np.linalg.eig(matrix.T@matrix)
            eigenvalues[iatom,jatom] = -np.sqrt(eigenvalues[iatom,jatom])

    # Get bond force constants
    bond_kAB = []
    for bond in bond_list:
        iatom=bond[0]; jatom=bond[1]
        AB = geometry[bond[1]] - geometry[bond[0]]
        AB /= np.linalg.norm(AB)
        k_AB = 0.0
        for i in range(3):
            k_AB -= eigenvalues[bond[0],bond[1],i]*np.abs(np.dot(AB,eigenvectors[bond[0],bond[1],:,i]))
            k_AB -= eigenvalues[bond[1],bond[0],i]*np.abs(np.dot(-AB,eigenvectors[bond[1],bond[0],:,i]))
        k_AB *= 0.5
        bond_kAB.append(k_AB)

    # Get angle force constants
    angle_kAB = []
    for i,triad1 in enumerate(angle_list):
        alist = []
        clist = []
        for j,triad2 in enumerate(angle_list):
            if i == j: continue
            if triad1[1] == triad2[1]:
                if (triad1[0] == triad2[0]) or (triad1[0] == triad2[2]):
                    alist.append(triad2)
                if (triad1[2] == triad2[0]) or (triad1[2] == triad2[2]):
                    clist.append(triad2)

        if len(alist) == 0:
            ascaling = 1.0
        else:
            ascaling = 0.0
            AB = geometry[triad1[1]] - geometry[triad1[0]]
            CB = geometry[triad1[1]] - geometry[triad1[2]]
            AB /= np.linalg.norm(AB)
            CB /= np.linalg.norm(CB)
            u_N = np.cross(CB, AB)
            u_N /= np.linalg.norm(u_N)
            PA1 = np.cross(u_N, AB)
            PA1 /= np.linalg.norm(PA1)
            for triad2 in alist:
                CB = copy.copy(geometry[triad1[2]])
                CB -= geometry[triad2[2]] if triad1[0] == triad2[0] else geometry[triad2[0]]
                CB /= np.linalg.norm(CB)
                u_N = np.cross(CB, AB)
                u_N /= np.linalg.norm(u_N)
                PA = np.cross(u_N, AB)
                PA /= np.linalg.norm(PA)
                ascaling += np.dot(PA1,PA)**2
            ascaling /= len(alist)
            ascaling += 1.0

        if len(clist) == 0:
            cscaling = 1.0
        else:
            cscaling = 0.0
            AB = geometry[triad1[1]] - geometry[triad1[0]]
            CB = geometry[triad1[1]] - geometry[triad1[2]]
            AB /= np.linalg.norm(AB)
            CB /= np.linalg.norm(CB)
            u_N = np.cross(AB, CB) 
            u_N /= np.linalg.norm(u_N)
            PA2 = np.cross(u_N, CB)
            PA2 /= np.linalg.norm(PA2)
            for triad2 in alist:
                AB = copy.copy(geometry[triad1[1]])
                AB -= geometry[triad2[2]] if triad1[2] == triad2[0] else geometry[triad2[0]]
                AB /= np.linalg.norm(AB)
                u_N = np.cross(AB, CB)
                u_N /= np.linalg.norm(u_N)
                PA = np.cross(u_N, CB)
                PA /= np.linalg.norm(PA)
                cscaling += np.dot(PA2,PA)**2
            cscaling /= len(clist)
            cscaling += 1.0
        
        k_AB = 0.0
        for ibond,bond in enumerate(bond_list):
            if (triad1[0] in bond) and (triad1[1] in bond):
                k_PA1 = 0.0
                for i in range(3):
                    k_PA1 -= eigenvalues[bond[0],bond[1],i]*np.abs(np.dot(PA1,eigenvectors[bond[0],bond[1],:,i]))
                    k_PA1 -= eigenvalues[bond[1],bond[0],i]*np.abs(np.dot(-PA1,eigenvectors[bond[1],bond[0],:,i]))
                k_PA1 *= 0.5
                k_AB += ascaling/(k_PA1 * bond_lengths[triad1[0],triad1[1]]**2)
            elif (triad1[2] in bond) and (triad1[1] in bond):
                k_PA2 = 0.0
                for i in range(3):
                    k_PA2 -= eigenvalues[bond[0],bond[1],i]*np.abs(np.dot(PA2,eigenvectors[bond[0],bond[1],:,i]))
                    k_PA2 -= eigenvalues[bond[1],bond[0],i]*np.abs(np.dot(-PA2,eigenvectors[bond[1],bond[0],:,i]))
                k_PA2 *= 0.5
                k_AB += cscaling/(k_PA2 * bond_lengths[triad1[0],triad1[1]]**2)

        angle_kAB.append(1.0/k_AB)

    return np.array(bond_kAB)*100.0, np.array(angle_kAB), bond_list, angle_list



def bond_angle_list(elements,geometry):
    _natoms = len(elements)
    _covrad = [ covrad[element] for element in elements]

    bonds = []
    # find bonded pairs
    for iatom in range(_natoms):
        for jatom in range(iatom+1,_natoms):
            dist = np.linalg.norm(geometry[iatom]-geometry[jatom])
            if dist < 1.2*(_covrad[iatom] + _covrad[jatom]):
                bonds.append([iatom,jatom])
    bonds = sorted(bonds, key=lambda x:(int(x[0]),int(x[1])))
    bonds = np.array(bonds,dtype=int)

    # find bonded triplets
    angles = []
    for ibond in range(len(bonds)):
        for jbond in range(ibond+1,len(bonds)):
            if bonds[ibond,0] == bonds[jbond,0]:
                _first  = bonds[ibond,1]
                second = bonds[ibond,0]
                _third  = bonds[jbond,1]
            elif bonds[ibond,0] == bonds[jbond,1]:
                _first  = bonds[ibond,1]
                second = bonds[ibond,0]
                _third  = bonds[jbond,0]
            elif bonds[ibond,1] == bonds[jbond,0]:
                _first  = bonds[ibond,0]
                second = bonds[ibond,1]
                _third  = bonds[jbond,1]
            elif bonds[ibond,1] == bonds[jbond,1]:
                _first  = bonds[ibond,0]
                second = bonds[ibond,1]
                _third  = bonds[jbond,0]
            else:
                continue

            first = min(_first,_third)
            third = max(_first,_third)

            angles.append([first,second,third])

    angles = sorted(angles, key=lambda x:(int(x[1]),int(x[0]),int(x[2])))
    return bonds,np.array(angles,dtype=int)




