from ptmpsi.nwchem.templates import *
import copy
import numpy as np

def readoptim(filename):
    geometry = ""
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
            found = (line.find("Optimization converged") > -1 or
                     line.find("Failed to converge in ") > -1)
            if found: break
            line = fh.readline()

        if not found:
            print(" Error: could not find last geometry")
            quit()

        found = False
        line = fh.readline()
        while line:
            found = line.find("Output coordinates") > -1
            if found: break
            line = fh.readline()

        if not found:
            print(" Error: could not find last geometry")
            quit()

        fh.readline()
        fh.readline()
        fh.readline()
        line = fh.readline()
        while line:
            if not line.strip(): break
            split = line.split()
            _natoms = int(split[0])
            geometry += coordinates.format(split[1],*[float(x) for x in split[3:]])
            line = fh.readline()

        if natoms != _natoms:
            print(" Error: number of atoms mismatch")
            quit()
    return natoms,charge,mult,geometry


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
