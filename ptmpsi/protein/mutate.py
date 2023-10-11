import numpy as np
from copy import deepcopy as copy
from ptmpsi.math import alignres, rotate_chi1, rotate_chi2, find_clashes_residue, nerf, rotmatvec
from ptmpsi.residues import Residue, resdict, ptmdict, ptm2nonstandard
from ptmpsi.residues.ptms import doptm, check_ptm, get_ptm_name, add_hydrogens
from ptmpsi.residues.template import Template
from ptmpsi.exceptions import MyDockingError
from ptmpsi.constants import nhbond, amidebond
from ptmpsi.protein.tools import get_residue, get_template

_CYSPTMS = ["carbamoylation", "sulfhydration","sulfenylation","sulfinylation","sulfonylation","nitrosylation","glutathionylation", "cysteinylation"]

def point_mutation(protein,original,new):
    """
    Perform a point mutation.
    """
    # Get Residue and Template instances
    _original = get_residue(protein,original) 
    _new = get_template(new)

    # Get new coordinates
    newcoords = alignres(_original,_new)

    # Update residue fields
    oldatoms = copy(_original.natoms)
    _original.coordinates = copy(newcoords)
    _original.names = copy(_new.elements[:,1])
    _original.elements = copy(_new.elements[:,0])
    _original.backbone = copy(_new.backbone)
    _original.chi1 = copy(_new.chi1)
    _original.chi2 = copy(_new.chi2)
    _original.cattach = copy(_new.cattach)
    _original.nattach = copy(_new.nattach)
    _original.natoms = len(newcoords)
    _original.name = _new.name

    # Fix amide hydrogen position
    if _original.resid > 1:
        for chain in protein.chains:
            if chain.name == _original.chain:
                previous = chain.residues[_original.resid-2]
                h = nerf(previous.find_coord("O"),
                         previous.find_coord("C"),
                         _original.find_coord("N"),
                         nhbond, 120, 180)
                _original.coordinates[_original.find("H")] = copy(h)
        newcoords = copy(_original.coordinates)

    # Update protein
    if oldatoms != _original.natoms:
        protein.update()

    # Find possible clashes
    found = True
    nclashes = find_clashes_residue(_original,[protein],printing=False)

    # Scan chi1 and chi2 dihedrals for a better rotamer
    if nclashes > 0:
        found, angle1, angle2, minclashes = scan_chi1_chi2(protein,_original,nclashes,_new.chi1,_new.chi2)

    # Get final rotamer
    if found:
        print("\n\t Found rotamer with no clashes!")
    else:
        print("\n\t\t Warning: All rotamers had clashes!")
        print("\t\t          using rotamer with {} clashes".format(minclashes))
        _original.coordinates = copy(newcoords)
        if angle1 > 0: rotate_chi1(_original,_new.chi1,angle1)
        if angle2 > 0: rotate_chi2(_original,_new.chi2,angle2)



            
    return


def post_translational_modification(protein,original,ptm):
    _ptm = ptm.lower()
    _original = get_residue(protein,original)

    # Check if residue is either C- or N-terminus
    for chain in protein.chains:
        if chain.name == _original.chain:
            nterminus = _original == chain.residues[0]
            cterminus = _original == chain.residues[-1]
            # Check for N-terminal acetylation case
            if nterminus and (ptm == "alpha-acetylation"):
                protein.prepend(_original.chain,"ACE")
                return

    # Get radical to be attached
    _radical = ptmdict.get(_ptm,0)
    if _radical == 0:
        raise MyDockingError("Post-translational modification '{}' is not coded".format(_ptm))

    # Special case for Cystein PTMs
    if _original.name in ["CYS", "CYS", "CYM"]:
        new = ptm2nonstandard.get(_ptm,None)
        print(ptm,new)
        if new is not None:
            point_mutation(protein,_original,new)
            return
    elif _ptm in _CYSPTMS:
        raise MyDockingError("Post-translational modification '{}' is only coded for CYS-type residues".format(_ptm))

    # TODO:
    if _radical is None:
        raise MyDockingError("TODO")

    # Check if PTM is compatible with residue and get geometrical parameters
    bond, angle, dihedral = check_ptm(_original.name,_ptm)

    # Perform PTM
    doptm(_original,_radical,bond,angle,dihedral)
    if _ptm == 'dimethylation': 
        if _original.name in ["LYS","LYN"]:
            doptm(_original,_radical,bond,109,150)
        else:
            doptm(_original,_radical,bond,120,0)
    elif _ptm == 'trimethylation':
        doptm(_original,_radical,bond,109,150)
        doptm(_original,_radical,bond,109,30)
    elif _ptm == "symmetric dimethylation":
        doptm(_original,_radical,bond,120,180,True)
    elif _ptm == "asymmetric dimethylation":
        doptm(_original,_radical,bond,120,0,True)

    # Add missing hydrogens
    if _original.name  in ["ARG","LYS","LYN"]:
        add_hydrogens(_original,_ptm)

    newcoords = copy(_original.coordinates)

    # Update protein
    protein.update()

    # Find possible clashes
    found = True
    nclashes = find_clashes_residue(_original,[protein])

    # Scan chi1 and chi2 dihedrals for a better rotamer
    if nclashes > 0:
        _template = copy.deepcopy(resdict[_original.name])
        if _template.chi1 is None:
            chi1 = None
        else:
            chi1 = np.array([
                _original.find(_template.elements[_template.chi1[0],1]),
                _original.find(_template.elements[_template.chi1[1],1]),
                _original.find(_template.elements[_template.chi1[2],1]),
                _original.find(_template.elements[_template.chi1[3],1])
                ])

        if _template.chi2 is None:
            chi2 = None
        else:
            chi2 = np.array([
                _original.find(_template.elements[_template.chi2[0],1]),
                _original.find(_template.elements[_template.chi2[1],1]),
                _original.find(_template.elements[_template.chi2[2],1]),
                _original.find(_template.elements[_template.chi2[3],1])
                ])
        found, angle1, angle2, minclashes = scan_chi1_chi2(protein,_original,nclashes,chi1,chi2)

    # Get final rotamer
    if found:
        print("\n\t Found rotamer with no clashes!")
    else:
        print("\n\t\t Warning: All rotamers had clashes!")
        print("\t\t          using rotamer with {} clashes".format(minclashes))
        _original.coordinates = copy(newcoords)
        if angle1 > 0: rotate_chi1(_original,chi1,angle1)
        if angle2 > 0: rotate_chi2(_original,chi2,angle2)

    # Update residue name
    _original.name = get_ptm_name(_original.name, _ptm)

    return


def scan_chi1_chi2(protein,residue,nclashes,chi1,chi2):
    print("\n\t Current rotamer has {} possible clashes".format(nclashes))
    angle1 = 0; angle2 = 0
    minclashes = nclashes
    found = False

    # If chi1 info is None, no possible rotamers (GLY, PRO)
    if (chi1 is not None):
        rotamers1 = range(0,360,30)
        rotamers2 = [0] if chi2 is None else range(0,360,30)
        for irot in rotamers1:
            if found: break
            if irot > 0:
                internal = rotate_chi1(residue,chi1,30)
            for jrot in rotamers2:
                if (irot == 0) and (jrot == 0): continue
                if chi2 is not None:
                    internal = rotate_chi2(residue,chi2,30)
                if internal: continue
                nclashes = find_clashes_residue(residue,[protein],printing=False)
                if nclashes == 0:
                    found = True
                    return found, 0, 0, 0
                print("\t Current rotamer has {} possible clashes".format(nclashes))
                if nclashes < minclashes:
                    minclashes = nclashes
                    angle1 = irot
                    angle2 = jrot

    return found, angle1, angle2, minclashes
