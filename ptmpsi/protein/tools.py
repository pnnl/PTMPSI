import copy
import numpy as np
from itertools import combinations as get_combinations
from ptmpsi.protein import Residue
from ptmpsi.residues import resdict, three2one, one2three, ptm2nonstandard, ptmdict
from ptmpsi.residues.template import Template
from ptmpsi.residues.ptms import check_ptm
from ptmpsi.exceptions import MyDockingError
from ptmpsi.gromacs import generate as generate_gromacs
from ptmpsi.gromacs.utils import amber_to_gromacs_names

_CYSPTMS = ["carbamoylation", "sulfhydration","sulfenylation","sulfinylation","sulfonylation","nitrosylation","glutathionylation", "cysteinylation"]

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
        else:
            lower = 1
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
            if len(residue) == 1:
                _residue = copy.deepcopy(resdict[one2three[residue.upper()]])
            else:
                _residue = copy.deepcopy(resdict[residue.upper()])
        except:
            raise MyDockingError("There is no template for '{}'".format(residue))
    else:
        raise MyDockingError("New residue should be a string or a Template object")
    return _residue

def findresidue(protein, name):
    try:
        residues = [ f"{chain.name}:{residue.name}{residue.resid}" for chain in protein.chains for residue in chain.residues if residue.name == name ]
    except:
        residues = [ f"{protein.name}:{residue.name}{residue.resid}" for residue in protein.residues if residue.name == name ]
    return residues


def check_res_str(residue):
    _chain = _resnum = _resname = None

    # Check if a generic residue name was given
    oneletter = three2one.get(residue, None)
    if oneletter is not None:
        return _chain, residue, _resnum

    threeletter = one2three.get(residue, None)
    if threeletter is not None:
        return _chain, threeletter, _resnum

    # Check if chain was given 
    _chain = None
    if ":" in residue[1]:
        _chain = residue[0]
        _residue = residue[2:].upper()
    else:
        _residue = residue.upper()

    oneletter = three2one.get(_residue, None)
    if oneletter is not None:
        return _chain, _residue, _resnum

    threeletter = one2three.get(_residue, None)
    if threeletter is not None:
        return _chain, threeletter, _resnum

    # 3-letter aminoacid code
    if _residue[1].isalpha():
        if len(_residue) < 4:
            raise KeyError("'{}' is not a 3-letter code + residue number".format(residue))
        _lower = 3
    elif _residue[0].isalpha():
        _lower = 1
    else:
        _lower = 0

    # Check that the rest of the characters are digits
    for digit in _residue[_lower:]:
        assert digit.isdigit(), f"'{_residue[lower:]}' is not a number"

    _resnum = int(_residue[_lower:])
    if _lower == 1:
        _resname = one2three(_residue[0])
    elif _lower == 3:
        _resname = _residue[:_lower]
    
    return _chain, _resname, _resnum


def ptm_combination(protein, ptms, ntuple=-1, exclude=None):
    assert isinstance(ptms, list), "ptms must be a List of tuples (RES, PTM)"
    assert len(ptms) > 0, "ptms must have size > 0"

    # Check that all provided PTMs are valid and compute total number of combinations
    ptmlist = []
    for res, ptm in ptms:
        assert isinstance(res, str), f"Residue name must be a string of the type\n\tC:RRRNNN\n\tRRRNNN\n\tNNN"
        assert isinstance(ptm, str), f"ptm must be a string"

        _ptm = ptm.lower()
        _radical = ptmdict.get(_ptm, 0)
        assert _radical != 0, f"Post-translational modification '{_ptm}' is not available"

        _chain, _resname, _resnum = check_res_str(res)
        if _chain is None and _resnum is None:
            residues = findresidue(protein, _resname)
        elif _resnum is None:
            chain = [ c for c in protein.chains if c.name == _chain ]
            residues = findresidue(chain[0], _resname)
        elif _chain is None and _resname is None:
            if (len(protein.chains) == 1):
                residues = [ f"{protein.chains[0].name}:{protein.chains[0].residues[_resnum-1].name}{_resnum}" ]
            else:
                residues = [ f"{c.name}:{r.name}{_resnum}" for c in protein.chains for r in c.residues ]
            _resname = residues[0][2:4]
            for r in residues:
                assert  r[2:4] == _resname, f"Ambiguous residue specification '{res}'" 
        elif _resname is None:
            chain = [ c for c in protein.chains if c.name == _chain ]
            _resname = chain[0].residues[_resnum-1].name
            residues = [ f"{_chain}:{_resname}{_resnum}" ]
        elif _chain is not None:
            chain = [ c for c in protein.chains if c.name == _chain ]
            assert _resname == chain[0].residues[_resnum-1].name, f"Residue '{res}' does not exist: {chain[0].name}:{chain[0].residues[_resnum-1].name}{_resnum}"
            residues = [ f"{chain[0].name}:{_resname}{_resnum}" ]
        else:
            if len(protein.chains) == 1:
                residues = [ f"{protein.chains[0].name}:{_resname}{_resnum}" ]
            else:
                residues = [ f"{c.name}:{r.name}{_resnum}" for c in protein.chains for r in c.residues ]
                for r in residues:
                    assert r[2:4] == _resname, f"Ambiguous residue specification '{res}'" 

        _dummy = None
        if _resname in ["CYS", "CYM", "CYX"]:
            _dummy = ptm2nonstandard.get(_ptm, None)
        else:
            assert _ptm not in _CYSPTMS, f"Post-translational modification '{_ptm}' is only coded for CYS-type residues"

        if _dummy is None:
            _dummy, _dummy, _dummy = check_ptm(_resname, _ptm)

        # Up to this point we know that the RES: PTM combination is valid
        for _residue in residues:
            ptmlist.append( (_residue, _ptm) )

    # Remove duplicates
    ptmlist = list(dict.fromkeys(ptmlist))

    # Check for exclusions
    if exclude is not None:
        assert isinstance(exclude, list), "exclude must be a List of residues with format C:RRRNNN"
        ptmlist = [ x for x in ptmlist if x[0] not in exclude ]
        
    # Get number of cases
    levels = len(ptmlist) if ntuple == -1 else ntuple
    assert levels > 0, "ntuple should be greater than 0"

    nperm = []
    combinations = []
    for level in range(levels):
        _combinations = [ *get_combinations(ptmlist, level+1)]
        _valid = []
        for _combination in _combinations:
            _residues = [ x[0] for x in _combination ]
            _test = set(_residues)
            if len(_test) < len(_residues): continue
            _valid.append(_combination)
        combinations.append( _valid.copy() )
        nperm.append( len(combinations[level]) )
        print(f"\tNumber of {level+1}-tuple combinations: {nperm[level]}")
    print(f"Total number of combinations: { np.sum(np.array(nperm)) }")

    return combinations

