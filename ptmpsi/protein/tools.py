from ptmpsi.protein import Residue
from ptmpsi.residues import resdict
from ptmpsi.residues.template import Template
from ptmpsi.exceptions import MyDockingError

def get_residue(protein,residue):
    # If residue is a residue instance, nothing to be done
    if isinstance(residue,Residue):
        _residue = residue

    # If residue is a string
    elif isinstance(residue,str):
        # Check if chain was given
        chain = None
        if ":" in residue[1]:
            chain = residue[0]
            residue = residue[2:].upper()
        # 3-letter aminoacid code
        if residue[1].isalpha():
            if len(residue) < 4:
                raise KeyError("'{}' is not a 3-letter code + residue number".format(residue))
            lower = 3
        # Check that the rest of the characters are digits
        for digit in residue[lower:]:
            if not digit.isdigit():
                raise KeyError("'{}' is not a residue number".format(residue[lower:]))
        resnum = int(residue[lower:])
        # Get the Residue instance
        if chain is not None:
            _residue = None
            for _chain in protein.chains:
                if _chain.name == chain:
                    _residue = _chain.residues[resnum-1]
                    break
            if _residue is None: raise MyDockingError("There is no chain '{}' in protein".format(chain))
        else:
            _residue = protein.chains[0].residues[resnum-1]
        # Check that the name and number correspond
        if lower == 3:
            if _residue.name != residue[:3]:
                raise KeyError("Residue '{}' could not be found".format(residue))

    # If residue is an integer, assume first chain
    elif isinstance(residue,int):
        resnum = int(residue)
        _residue = protein.chains[0].residues[resnum-1]
    else:
        raise MyDockingError("New residue should be a string or a Residue object")

    return _residue


def get_template(residue):
    # Check that the
    if isinstance(residue,Template):
        _residue = residue
    elif isinstance(residue,str):
        try:
            _residue = resdict[residue.upper()]
        except:
            raise MyDockingError("There is no template for '{}'".format(residue))
    else:
        raise MyDockingError("New residue should be a string or a Template object")
    return _residue
