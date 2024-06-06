import numpy as np
from ptmpsi.polymers.tools import UnitCell, gen_pdb, gen_xyz

alphabet = {"CA":6,"CB":5, "CG":4, "CD":3, "CE":2, "N":1, "C":0}
dict2 = {"N": "CE", "CE":"CD", "CD":"CG", "CG":"CB", "CB":"CA", "CA":"C", "C":"O", "O":"N"}

models = {}
models["Malta1979"] = UnitCell(
  9.71, 8.19, 17.24, 90.0, 90.0, 115.0, 
   np.array([
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
   ]),
   ref="EPJ 15, 765 (1979)")

models["Holmes"] = UnitCell(
  9.56, 8.01, 17.24, 90.0, 90.0, 112.5,
np.array([
["N",  0.037,  0.006,  0.285],
["C", -0.040, -0.007,  0.219],
["O", -0.160, -0.027,  0.219],
["C",  0.044,  0.007,  0.146],
["C", -0.040, -0.007,  0.073],
["C",  0.044,  0.007,  0.000],
["C", -0.040, -0.007, -0.073],
["C",  0.044,  0.007, -0.146],
["N",  0.537,  0.006,  0.215],
["C",  0.460, -0.007,  0.281],
["O",  0.340, -0.027,  0.281],
["C",  0.544,  0.007,  0.354],
["C",  0.460, -0.007,  0.427],
["C",  0.544,  0.007,  0.500],
["C",  0.460, -0.007,  0.573],
["C",  0.544,  0.007,  0.646],
["N",  0.037,  0.506,  0.500],
["C", -0.040,  0.493,  0.433],
["O", -0.160,  0.473,  0.433],
["C",  0.044,  0.507,  0.360],
["C", -0.040,  0.493,  0.287],
["C",  0.044,  0.507,  0.214],
["C", -0.040,  0.493,  0.141],
["C",  0.044,  0.507,  0.068],
["N",  0.537,  0.506,  0.429],
["C",  0.460,  0.493,  0.495],
["O",  0.340,  0.473,  0.495],
["C",  0.544,  0.507,  0.568],
["C",  0.460,  0.493,  0.641],
["C",  0.544,  0.507,  0.714],
["C",  0.460,  0.493,  0.787],
["C",  0.544,  0.507,  0.860]
]),
 ref="Holmes"
)


models["Arimoto"] = UnitCell(
  9.33, 4.78, 16.88, 90.0, 90.0, 121.0,
np.array([
["N",  0.306,  0.114,  0.022],
["C",  0.290,  0.116,  0.169],
["O",  0.236, -0.394, -0.016],
["C",  0.390,  0.108,  0.097],
["C",  0.236, -0.138, -0.029],
["C",  0.153, -0.105, -0.103],
["C",  0.236, -0.109, -0.181],
["C",  0.134, -0.105, -0.251],

["N",  0.537,  0.006,  0.215],
["C",  0.460, -0.007,  0.281],
["O",  0.340, -0.027,  0.281],
["C",  0.544,  0.007,  0.354],
["C",  0.460, -0.007,  0.427],
["C",  0.544,  0.007,  0.500],
["C",  0.460, -0.007,  0.573],
["C",  0.544,  0.007,  0.646],
["N",  0.037,  0.506,  0.500],
["C", -0.040,  0.493,  0.433],
["O", -0.160,  0.473,  0.433],
["C",  0.044,  0.507,  0.360],
["C", -0.040,  0.493,  0.287],
["C",  0.044,  0.507,  0.214],
["C", -0.040,  0.493,  0.141],
["C",  0.044,  0.507,  0.068],
["N",  0.537,  0.506,  0.429],
["C",  0.460,  0.493,  0.495],
["O",  0.340,  0.473,  0.495],
["C",  0.544,  0.507,  0.568],
["C",  0.460,  0.493,  0.641],
["C",  0.544,  0.507,  0.714],
["C",  0.460,  0.493,  0.787],
["C",  0.544,  0.507,  0.860]
]),
 ref="Holmes"
)


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

def gen_surface(rangea=1, rangeb=1, rangec=1, model="Malta1979"):
    cell = models[model]
    xyzcoords, labels = gen_xyz(cell, rangea, rangeb, rangec)
    gen_pdb(xyzcoords, labels, alphabet, dict2, rename, "nylon6", f"{rangea}{rangeb}{rangec}")
    return

