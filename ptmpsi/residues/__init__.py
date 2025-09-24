import numpy as np
from . import aminoacids
from . import nonstandard
from . import ptms
from ptmpsi.exceptions import MyDockingError

class Atom:
    def __init__(self, name=None, coord=None, bfactor=0.0, occupancy=1.0, altloc=None, element=None):
        if coord is None:
            raise MyDockingError(f"Invalid coord specification '{coord}'. Failed to initialize Atom instance")
        elif isinstance(coord, (list, np.ndarray)):
            if len(coord) != 3:
                raise MyDockingError(f"Invalid coord specification '{coord}'. Failed to initialize Atom instance")
        self.name = name if name is not None else "DU"
        self.element = element if element is not None else "H"
        self.bfactor = bfactor
        self.occupancy = occupancy
        self.altloc = altloc
        self.coords = coords
            

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
        self.chi1  = None
        self.chi2  = None
        self.nattach = np.empty(3,dtype=float)
        self.cattach = np.empty(3,dtype=float)

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

    def find_coord(self, atom):
        return self.coordinates[self.find(atom)]

    def add(self, atom):
        if not isinstance(atom, Atom):
            raise MyDockingError("atom argument must be of Atom class")
        self.natoms += 1
        names = np.empty(self.natoms, dtype='U4')
        elements = np.empty(self.natoms, dtype='U4')
        coordinates = np.empty((self.natoms, 3), dtype=float)
        for i in range(self.natoms-1):
            names[i] = self.names[i]
            elements[i] = self.elements[i]
            coordinates[i] = self.coordinates[i]
        names[-1] = atom.name
        elements[-1] = atom.element
        coordinates[-1] = atom.coords
        self.names = copy.deepcopy(names)
        self.elements = copy.deepcopy(elements)
        self.coordinates = copy.deepcopy(coordinates)



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
        "HIS": "H",
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
        "reduction": None,
        "oxidation": None,
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
        "cysteinylation": "IYY",
        "reduction": "CYS",
        "oxidation": "CYX",
        }
