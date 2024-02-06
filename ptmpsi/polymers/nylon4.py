import numpy as np
from ptmpsi.polymers.tools import UnitCell, gen_pdb, gen_xyz

alphabet = {"CA":4,"CB":3, "CG":2, "N":1, "C":0}
dict2 = {"N": "CG", "CG":"CB", "CB":"CA", "CA":"C", "C":"O", "O":"N"}

models = {}
models["Fredericks"] = UnitCell(
  9.29, 7.27, 12.24, 90.0, 90.0, 114.5, 
   np.array([
   ["C",-0.044000,0.507000,0.300000],
   ["C", 0.040000,0.493000,0.402000],
   ["C",-0.044000,0.507000,0.504000],
   ["N",-0.037000,0.506000,0.699000],
   ["C", 0.040000,0.493000,0.607000],
   ["O", 0.160000,0.473000,0.607000],
   ["C", 0.044000,1-0.507000,0.800000],
   ["C",-0.040000,1-0.493000,0.902000],
   ["C", 0.044000,1-0.507000,1.004000],
   ["N", 0.037000,1-0.506000,1.199000],
   ["C",-0.040000,1-0.493000,1.107000],
   ["O",-0.160000,1-0.473000,1.107000],
#    
   ["C",-0.044000,0.007000,0.500000],
   ["C", 0.040000,0.007000,0.602000],
   ["C",-0.044000,0.007000,0.704000],
   ["N",-0.037000,0.006000,0.899000],
   ["C", 0.040000,0.007000,0.807000],
   ["O", 0.160000,0.027000,0.807000],
   ["C", 0.044000,-0.007000,1.000000],
   ["C",-0.040000,-0.007000,1.102000],
   ["C", 0.044000,-0.007000,1.204000],
   ["N", 0.037000,-0.006000,1.399000],
   ["C",-0.040000,-0.007000,1.307000],
   ["O",-0.160000,-0.027000,1.307000],

#    
   ["C", 0.460000,0.007000,0.500000],
   ["C", 0.544000,0.007000,0.602000],
   ["C", 0.460000,0.007000,0.704000],
   ["N", 0.537000,0.006000,0.807000],
   ["C", 0.460000,0.007000,0.899000],
   ["O", 0.340000,0.027000,0.899000],
   ["C", 1-0.460000,0.007000,1.000000],
   ["C", 1-0.544000,0.007000,1.102000],
   ["C", 1-0.460000,0.007000,1.204000],
   ["N", 1-0.537000,0.006000,1.307000],
   ["C", 1-0.460000,0.007000,1.399000],
   ["O", 1-0.340000,0.027000,1.399000],
#    
   ["C", 0.540000,0.507000,0.300000],
   ["C", 0.456000,0.493000,0.402000],
   ["C", 0.540000,0.507000,0.504000],
   ["N", 0.463000,0.494000,0.607000],
   ["C", 0.540000,0.507000,0.699000],
   ["O", 0.660000,0.527000,0.699000],
   ["C", 1-0.540000,1-0.507000,0.800000],
   ["C", 1-0.456000,1-0.493000,0.902000],
   ["C", 1-0.540000,1-0.507000,1.004000],
   ["N", 1-0.463000,1-0.494000,1.107000],
   ["C", 1-0.540000,1-0.507000,1.199000],
   ["O", 1-0.660000,1-0.527000,1.199000],
]),
   ref="EPJ 15, 765 (1979)")

def rename(coords,bonds,natoms):
    for iatom in range(natoms):
        if coords[iatom,0] == 'C':
            obond = False
            nbond = False
            for jatom in range(natoms):
                if bonds[iatom,jatom]:
                    if coords[jatom,0] == 'O': obond = True
                    if coords[jatom,0] == 'N': nbond = True
            if nbond and obond:
                coords[iatom,0] = 'CO'
            elif nbond:
                coords[iatom,0] = 'CG'
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
        if coords[iatom,0] == 'CO':
            coords[iatom,0] = 'C'
            


def gen_surface(rangea=1, rangeb=1, rangec=1, model="Fredericks"):
    cell = models[model]
    xyzcoords, labels = gen_xyz(cell, rangea, rangeb, rangec, False)
    gen_pdb(xyzcoords, labels, alphabet, dict2, rename, "nylon4", f"{rangea}{rangeb}{rangec}")
    return

