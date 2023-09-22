#!/usr/bin/env python

import numpy as np
import copy

ang2bohr = 1.8897259886
eh2kjmol = 2625.5002
nm2bohr  = ang2bohr*10.0

covrad = {{
        "H": 0.31,
        "Li": 1.28,
        "Be": 0.96,
        "B": 0.84,
        "C": 0.76,
        "N": 0.71,
        "O": 0.66,
        "F": 0.57,
        "Na": 1.66,
        "Mg": 1.41,
        "Al": 1.21,
        "Si": 1.11,
        "P": 1.07,
        "S": 1.05,
        "Cl": 1.02,
        "K": 2.03,
        "Ca": 1.76,
        "V": 1.53,
        "Cr": 1.39,
        "Mn": 1.61,
        "Fe": 1.52,
        "Co": 1.50,
        "Cu": 1.32,
        "Zn": 1.22,
        "Se": 1.20,
        }}

mass = {{
        "H": 1.00797,
        "Li": 7.0,
        "Be": 9.012183,
        "B": 10.81,
        "C": 12.011,
        "N": 14.0067,
        "O": 15.9994,
        "F": 18.99840316,
        "Na": 22.9897693,
        "Mg": 24.305,
        "Al": 26.981538,
        "Si": 28.085, 
        "P": 30.97376200,
        "S": 32.07,
        "Cl": 35.45,
        "K": 39.0983,
        "Ca": 40.08,
        "V": 50.9415,
        "Cr": 51.996,
        "Mn": 54.93804,
        "Fe": 55.84,
        "Co": 58.93319,
        "Cu": 63.55,
        "Zn": 65.4,
        "Se": 78.97,
        }}

def read_xyz(filename):
    with open(filename,"r") as xyz:
        _natoms = int(xyz.readline())
        xyz.readline()
        geometry = np.zeros((_natoms,3))
        for i in range(_natoms):
            line = xyz.readline().split()[1:]
            geometry[i] = np.array(line, dtype=float)
    return geometry


def main():

    elements = {elements}
    for i in range(len(elements)):
        elements[i] = elements[i].capitalize()

    names    = {names}
    kabs,kthetas,alists,blists = [],[],[],[]

    for idx in range({nconf}):
        tail     = f"conf{{str(idx)}}"
        filename = f"{{tail}}.xyz"
        geometry = read_xyz(filename)
        kab,ktheta,blist,alist = bondedks(tail,elements,geometry)
        kabs.append(kab)
        kthetas.append(ktheta)
        alists.append(alist)
        blists.append(blist)

    kab = np.zeros(kabs[0].shape)
    for _kab in kabs:
        kab += _kab/len(kabs)

    ktheta = np.zeros(kthetas[0].shape)
    for _ktheta in kthetas:
        ktheta += _ktheta/len(kthetas)

    print(" Bond force constants [kJ/mol/nm^2]")
    print(" Atom1     Atom2       kAB")
    print(" ------------------------------")
    for i in range(len(blists[0])):
        print("{{: 3d}} {{:5s}}    {{: 3d}} {{:5s}}    {{:10.1f}}".format(
            blists[0][i,0],names[blists[0][i,0]],
            blists[0][i,1],names[blists[0][i,1]],kab[i]))

    print("")
    print(" Angle force constants [kJ/mol/rad^2]")
    print(" Atom1     Atom2     Atom3      kAB")
    print(" --------------------------------------")
    for i in range(len(alists[0])):
        print(" {{: 3d}} {{:5s}}    {{: 3d}} {{:5s}}    {{: 3d}} {{:5s}}    {{:8.3f}}".format(
            alists[0][i,0],names[alists[0][i,0]],
            alists[0][i,1],names[alists[0][i,1]],
            alists[0][i,2],names[alists[0][i,2]],ktheta[i]))
    print("")
    return

def bondedks(filename,elements,geometry):
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

        AB = geometry[triad1[1]] - geometry[triad1[0]]
        CB = geometry[triad1[1]] - geometry[triad1[2]]
        AB /= np.linalg.norm(AB)
        CB /= np.linalg.norm(CB)
        u_N = np.cross(CB, AB)
        u_N /= np.linalg.norm(u_N)
        PA1 = np.cross(u_N, AB)
        PA1 /= np.linalg.norm(PA1)
        u_N = np.cross(AB, CB)
        u_N /= np.linalg.norm(u_N)
        PA2 = np.cross(u_N, CB)
        PA2 /= np.linalg.norm(PA2)

        if len(alist) == 0:
            ascaling = 1.0
        else:
            ascaling = 0.0
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
            CB = geometry[triad1[1]] - geometry[triad1[2]]
            CB /= np.linalg.norm(CB)
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
                k_AB += cscaling/(k_PA2 * bond_lengths[triad1[2],triad1[1]]**2)

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
            if dist < 1.3*(_covrad[iatom] + _covrad[jatom]):
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


def readhess(filename):
    natoms = None
    with open(filename,"r") as fh:
        found = False
        line = fh.readline()
        while line:
            if natoms is None:
                if line.find("No. of atoms     :") > -1:
                    natoms = int(line.split()[4])
                    fh.readline(); fh.readline(); fh.readline()
                    line = fh.readline()
                    charge = int(line.split()[2])
                    line = fh.readline()
                    mult = int(line.split()[2])
                    line = fh.readline()
                    continue
            found = line.find("MASS-WEIGHTED PROJECTED HESSIAN") > -1
            if found: break
            line = fh.readline()

        if not found:
            print(" Error: Could not find hessian")
            quit()

        nxyz = 3*natoms
        hessian = np.zeros((nxyz,nxyz))
        jxyz = 0

        fh.readline()
        while True:
            fh.readline(); fh.readline()
            line = fh.readline().split()
            nvalues = len(line)
            fh.readline()
            for ixyz in range(jxyz,nxyz):
                line = fh.readline().split()
                ival = 0
                for kxyz in range(jxyz,jxyz+nvalues):
                    if kxyz > ixyz: break
                    ival += 1
                    hessian[ixyz,kxyz] = float(line[ival].replace("D","E"))
            jxyz += nvalues
            if jxyz == nxyz: break
    for ixyz in range(nxyz):
        for jxyz in range(ixyz+1,nxyz):
            hessian[ixyz,jxyz] = copy.copy(hessian[jxyz,ixyz])
    return hessian, natoms

if __name__ == "__main__":
    main()
