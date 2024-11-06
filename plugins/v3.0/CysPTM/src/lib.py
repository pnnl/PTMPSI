# vim: set expandtab shiftwidth=4 softtabstop=4:

# === UCSF ChimeraX Copyright ===
# Copyright 2016 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  For details see:
# http://www.rbvi.ucsf.edu/chimerax/docs/licensing.html
# This notice must be embedded in or attached to all copies,
# including partial copies, of the software or any revisions
# or derivations thereof.
# === UCSF ChimeraX Copyright ===

_dependent_cache = {}
_independent_cache = {}

from chimerax.rotamers import RotamerLibrary, RotamerParams, UnsupportedResTypeError, NoResidueRotamersError

class CysPTMRotamerLibrary(RotamerLibrary):

    @property
    def display_name(self):
        return "CysPTM"

    @property
    def description(self):
        return "Post-translationally modified Cys rotamer library"

    @property
    def citation(self):
        return """Environmental Molecular Sciences Laboratory (2022)."""

    _rotamer_res_names = set(["QCS", "SMC", "OCS", "SNC", "CSD", "CSO", "CSS", "XCN", "CGL"])

    @property
    def residue_names(self):
        return self._rotamer_res_names

    @property
    def res_name_descriptions(self):
        mapping = super().res_name_descriptions
        mapping.update({
            "QCS": "S-carbamoyl cysteine",
            "SMC": "S-methyl cysteine",
            "OCS": "Cysteine sulfonic acid",
            "SNC": "S-nitroso cysteine",
            "CSD": "Cysteine sulfinic acid",
            "CSO": "Cysteine sulfenic acid",
            "CSS": "Cysteine persulfide",
            "XCN": "S-cyano cysteine",
            "CGL": "S-glutathionyl cysteine",
        })
        return mapping

    @property
    def res_name_mapping(self):
        return { "carbamoyl": "QCS", "methyl": "SMC", "sulfonic": "OCS", "nitroso": "SNC", "sulfinic": "CSD", "sulfenic": "CSO", "cyano": "XCN", "glutathionyl": "CGL" }

    def rotamer_params(self, res_name, phi, psi, *, cis=False):
        if phi is None or psi is None:
            file_name = res_name
            archive = "independentCysPTMRotamers.zip"
            cache = _independent_cache
        else:
            from math import floor
            phi = floor((phi + 5) / 10.0) * 10
            psi = floor((psi + 5) / 10.0) * 10
            file_name = "%s%d%d" % (res_name, phi, psi)
            archive = "dependentCysPTMRotamers.zip"
            cache = _dependent_cache
        return self._get_params(res_name, file_name, cache, archive)

    @property
    def ptm_template_func(self):
        from .lib import templateResidue
        return templateResidue


def templateResidue(session,res):
    from chimerax.atomic import AtomicStructure
    
    m = AtomicStructure(session)
    residue = m.new_residue(res, 'A', 1)
    
    from .data import coordInfo
    crd = coordInfo[res]

    from chimerax.atomic.struct_edit import add_atom
    from numpy import array
    for element, name, pos in crd:
        add_atom(name, element, residue, array(pos, dtype=float))


    m.connect_structure()


    return residue
