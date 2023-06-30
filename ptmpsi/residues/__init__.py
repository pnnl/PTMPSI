import numpy as np
from . import aminoacids
from . import nonstandard
from . import ptms
from ptmpsi.exceptions import MyDockingError

class Residue:
    def __init__(self, resname, natoms):
        self.name = resname
        self.natoms = natoms
        self.names = np.empty(self.natoms,dtype='U4')
        self.elements = np.empty(self.natoms,dtype='U4')
        self.coordinates = np.empty((self.natoms,3),dtype=float)
        self.backbone = np.empty(3,dtype=int)
        self.resid = None
        self.chain = None

    def __eq__(self,other):
        if (self.resid == other.resid) and (self.chain == other.chain) and (self.name == other.name):
            return True
        else:
            return False

    def find(self,atom):
        pos = next((idx for idx,val in np.ndenumerate(self.names) if val==atom),None)
        if pos is None:
            raise MyDockingError("Atom '{}' was not found in Residue '{}:{}{}'".format(atom,self.chain,self.name,self.resid))
        return pos[0]

    def find_coord(self,atom):
        return self.coordinates[self.find(atom)]



three2one = {
        "ACE": "X",
        "ALA": "A",
        "ARG": "R",
        "ASH": "N",
        "ASN": "N",
        "ASP": "D",
        "CYM": "C",
        "CYS": "C",
        "CYX": "C",
        "GLH": "Q",
        "GLN": "Q",
        "GLU": "E",
        "GLY": "G",
        "HID": "H",
        "HIE": "H",
        "HIP": "H",
        "ILE": "I",
        "LEU": "L",
        "LYN": "K",
        "LYS": "K",
        "MET": "M",
        "NHE": "X",
        "NME": "X",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
        }


one2three = {
        "A": "ALA",
        "R": "ARG",
        "N": "ASN",
        "D": "ASP",
        "B": "ASX",
        "C": "CYS",
        "E": "GLU",
        "Q": "GLN",
        "Z": "GLX",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "L": "LEU",
        "K": "LYS",
        "M": "MET",
        "F": "PHE",
        "P": "PRO",
        "S": "SER",
        "T": "THR",
        "W": "TRP",
        "Y": "TYR",
        "V": "VAL",
        }

resdict = { "ACE": aminoacids.ACE,
            "ALA": aminoacids.ALA,
            "ARG": aminoacids.ARG,
            "ASH": aminoacids.ASH,
            "ASN": aminoacids.ASN,
            "ASP": aminoacids.ASP,
            "CYM": aminoacids.CYM,
            "CYS": aminoacids.CYS,
            "CYX": aminoacids.CYX,
            "GLH": aminoacids.GLH,
            "GLN": aminoacids.GLN,
            "GLU": aminoacids.GLU,
            "GLY": aminoacids.GLY,
            "HID": aminoacids.HID,
            "HIE": aminoacids.HIE,
            "HIP": aminoacids.HIP,
            "HIS": aminoacids.HIP,
            "ILE": aminoacids.ILE,
            "LEU": aminoacids.LEU,
            "LYN": aminoacids.LYN,
            "LYS": aminoacids.LYS,
            "MET": aminoacids.MET,
            "NHE": aminoacids.NHE,
            "NME": aminoacids.NME,
            "PHE": aminoacids.PHE,
            "PRO": aminoacids.PRO,
            "SER": aminoacids.SER,
            "THR": aminoacids.THR,
            "TRP": aminoacids.TRP,
            "TYR": aminoacids.TYR,
            "VAL": aminoacids.VAL,
#
            "QCS": nonstandard.QCS,
            "XCN": nonstandard.XCN,
            "SMC": nonstandard.SMC,
            "SNC": nonstandard.SNC,
            "CSD": nonstandard.CSD,
            "CSO": nonstandard.CSO,
            "OCS": nonstandard.OCS,
            "CSS": nonstandard.CSS,
            "CGL": nonstandard.CGL,
            "EAC": nonstandard.EAC,
            "ABA": nonstandard.ABA,
            "IYY": nonstandard.IYY,
        }

ptmdict = {
        "acetylation": ptms.acetylation,
        "citrullination": None,
        "cysteinylation": None,
        "glutathionylation": None,
        "glycosylation": None,
        "hydroxylation": None,
        "methylation": ptms.methylation,
        "myristoylation": None,
        "nitration": None,
        "nitrosylation": None,
        "palmitoylation": None,
        "phosphorylation": ptms.phosphorylation,
        "prenylation": None,
        "sulfhydration": None,
        "sulfenylation": None,
        "sulfinylation": None,
        "sulfonylation": None,
        "dimethylation": ptms.methylation,
        "trimethylation": ptms.methylation,
        "symmetric dimethylation": ptms.methylation,
        "asymmetric dimethylation": ptms.methylation,
        "cyanylation": None,
        "carbamoylation": None,
        }

ptm2nonstandard = {
        "glutathionylation": "CGL",
        "nitrosylation": "SNC",
        "sulfhydration": "CSS",
        "sulfenylation": "CSO",
        "sulfinylation": "CSD",
        "sulfonylation": "OCS",
        "methylation": "SMC",
        "carbamoylation": "QCS",
        "cyanylation": "XCN",
        "cysteinylation": "IYY"
        }
