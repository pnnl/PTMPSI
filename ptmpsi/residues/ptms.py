import numpy as np
from ptmpsi.constants import amidebond, pobond, nhbond
from ptmpsi.math import nerf, rotmatvec, get_torsion, rotmataxis

class PTM:
    def __init__(self):
        self.name = None
        self.chi1 = None
        self.chi2 = None
        self.names = None
        self.natoms = None
        self.attach = None
        self.elements = None
        self.coordinates = None


acetylation = PTM()
acetylation.name = "acetylation"
acetylation.coordinates = np.array([
 [2.000001, 1.000000, -1.346410E-06],
 [2.000001, 2.090000,  1.211769E-07],
 [1.486264, 2.453849,  0.889824    ],
 [1.486259, 2.453852, -0.889820    ],
 [3.427420, 2.640795, -2.981008E-06],
 [4.390580, 1.877406, -6.602501E-06]
])
acetylation.elements = np.array(["H", "C", "H", "H", "C", "O"])
acetylation.names = np.array(["HH31", "CH3", "HH32", "HH33", "CP", "OP"])
acetylation.natoms = len(acetylation.coordinates)
acetylation.attach = nerf(acetylation.coordinates[0],
                     acetylation.coordinates[1],
                     acetylation.coordinates[4],
                     amidebond,120,180)
acetylation.attach -= acetylation.coordinates[4]
acetylation.coordinates -= acetylation.coordinates[4]


phosphorylation = PTM()
phosphorylation.name = "phosphorylation"
phosphorylation.coordinates = np.array([
[ 0.000000, 0.000000,  0.000000],
[ 0.000000, 1.590000,  0.000000],
[ 1.499078, 2.119968,  0.000000],
[-0.749494, 2.120032, -1.298239]
])
phosphorylation.elements = np.array(["O","P","O","O"])
phosphorylation.names = np.array(["O1","P","O2","O3"])
phosphorylation.natoms = 4
phosphorylation.attach = np.array([-0.749522,2.120016,1.298230])
phosphorylation.attach -= phosphorylation.coordinates[1]
phosphorylation.coordinates -= phosphorylation.coordinates[1]


methylation = PTM()
methylation.name = "methylation"
methylation.coordinates = np.array([
 [2.000001, 1.000000, -1.346410E-06],
 [2.000001, 2.090000,  1.211769E-07],
 [1.486264, 2.453849,  0.889824    ],
 [1.486259, 2.453852, -0.889820    ],
    ])
methylation.elements = np.array(["H","C","H","H"])
methylation.names = np.array(["HH31","CH3","HH32","HH33"])
methylation.natoms = 4
methylation.attach = np.array([3.427420, 2.640795, -2.981008E-06])
methylation.attach -= methylation.coordinates[1]
methylation.coordinates -= methylation.coordinates[1]


hydrogens = {
        "NH1": ["NH11","NH12"],
        "NH2": ["NH21","NH22"],
        "OD2": ["HD2"],
        "SG": ["HG"],
        "ND1": ["HD1"],
        "NE2": ["HE2"],
        "NZ": ["HZ1","HZ2","HZ3"],
        "OG": ["HG"],
        "OG1": ["HG1"],
        "OH": ["HO"],
        }


ptmsite = {
        "ARG": ["NE","CZ","NH1"],
        "ASH": ["OD1","CG","OD2"],
        "ASP": ["OD1","CG","OD2"],
        "CYM": ["CA","CB","SG"],
        "CYS": ["CA","CB","SG"],
        "pros": ["NE2","CE1","ND1"],
        "tele": ["ND1","CE1","NE2"],
        "LYN": ["CD","CE","NZ"],
        "LYS": ["CD","CE","NZ"],
        "SER": ["CA","CB","OG"],
        "THR": ["CA","CB","OG1"],
        "TYR": ["CD1","CE1","OH"],
        }

def doptm(residue,radical,bond,angle,dihedral,argdimeth=False):


    # Ask for attachment site
    if residue.name in ["HIS","HIP","HID","HIE"]:
        print("\n\t Which histidine position should react?")
        print("\t\t 1. Pros (pi) position")
        print("\t\t 2. Tele (tau) position")
        position = input("\t Selection: ").lower()
        if position not in ["1","2","pros","tele"]:
            raise MyDockingError("Could not understand selection '{}'".format(position))
        site = ptmsite["pros" if position in ["1","pros"] else "tele"]
    elif argdimeth:
        site = ["NE","CZ","NH2"]
    else:
        site = ptmsite[residue.name]


    # Remove attached hydrogen atom (if any)
    mask = [*range(0,residue.natoms)]
    for hname in hydrogens[site[2]]:
        try:
            hpos = residue.find(hname)
            mask.pop(hpos)
        except:
            pass

    if len(mask) < residue.natoms:
        residue.elements = residue.elements[mask]
        residue.names = residue.names[mask]
        residue.coordinates = residue.coordinates[mask]
        residue.natoms = len(mask)

    # Define attachment point
    pos1 = residue.find_coord(site[2])
    pos2 = nerf(residue.find_coord(site[0]), residue.find_coord(site[1]),
                pos1,bond,angle,dihedral) - pos1


    # Align radical
    R = -rotmatvec(pos2,radical.attach)
    newcoords = np.dot(radical.coordinates,R.T) + pos1 + pos2

    # Attach radical
    residue.coordinates = np.vstack((residue.coordinates,newcoords))
    residue.names = np.hstack((residue.names,radical.names))
    residue.elements = np.hstack((residue.elements,radical.elements))
    residue.natoms += radical.natoms

    # Add missing hydrogens



    # Check that amide and esther bonds are planar
    if radical.name in ["acetylation","phosphorylation"]:
        site3 = "CP" if radical.name == "acetylation" else "P"
        site4 = "OP" if radical.name == "acetylation" else "O1"
        dihedral = get_torsion(residue.find_coord(site[1]),
                residue.find_coord(site[2]),residue.find_coord(site3),
                residue.find_coord(site4))
        if abs(dihedral) > 5:
            dihedral = -np.radians(dihedral)
            R = rotmataxis(residue.find_coord(site3)-residue.find_coord(site[2]),dihedral)
            mask = [ i for i,x in enumerate(residue.names) if x in radical.names ]
            temp = residue.coordinates[mask] - residue.find_coord(site[2])
            temp = np.dot(temp,R.T) + residue.find_coord(site[2])
            residue.coordinates[mask] = temp
            dihedral = get_torsion(residue.find_coord(site[1]),
                residue.find_coord(site[2]),residue.find_coord(site3),
                residue.find_coord(site4))


    return


def check_ptm(residue,ptm):
    if ptm == 'phosphorylation':
        # Check that the residue can be phosphorylated
        if residue not in ["SER","THR","TYR","ARG","HIS",
                           "HIP","HID","HIE","LYS","LYN", 
                           "ASH","ASP","CYS","CYM"]:
            raise MyDockingError("Don't know how to phosphorylate {} residue".format(residue))
        return pobond, 120, 180
    elif ptm == 'acetylation':
        if residue not in ["LYS","LYN"]:
            raise MyDockingError("Don't know how to acetylate {} residue".format(residue))
        return amidebond, 120, 180
    elif ptm == 'methylation':
        if residue not in ["GLU","GLH","LYS","LYN","ARG","HIS","HIP","HID","HIE"]:
            raise MyDockingError("Don't know how to methylate {} residue".format(residue))
        if residue in ["LYS","LYN"]:
            angle = 109
            dihedral = 270
        else:
            angle = 120
            dihedral = 180
        return amidebond, angle, dihedral
    elif ptm == 'dimethylation':
        if residue not in ["LYS","LYN","ARG"]:
            raise MyDockingError("Don't know how to dimethylate {} residue".format(residue))
        if residue == "ARG":
            angle = 120
            dihedral = 180
        else:
            angle = 109
            dihedral = 270
        return amidebond, angle, dihedral
    elif ptm == 'trimethylation':
        if residue not in ["LYS","LYN"]:
            raise MyDockingError("Don't know how to trimethylate {} residue".format(residue))
        return amidebond, 109, 270
    elif ptm == "symmetric dimethylation":
        if residue not in ["ARG"]:
            raise MyDockingError("Don't know how to symmetric dimethylate {} residue".format(residue))
        return amidebond, 120, 180
    elif ptm == "asymmetric dimethylation":
        if residue not in ["ARG"]:
            raise MyDockingError("Don't know how to asymmetric dimethylate {} residue".format(residue))
        return amidebond, 120, 180
    return


def get_ptm_name(original,ptm):
    # TODO
    return "PTM"


def add_hydrogens(residue,ptm):
    if (residue.name in ["LYS","LYN"]) and (ptm != "trimethylation"):
        if ptm in ["methylation","dimethylation"]:
            angle = 109
            dihedral = 30
        else:
            angle = 120
            dihedral = 0

        newh = nerf(residue.find_coord("CD"),residue.find_coord("CE"),
                        residue.find_coord("NZ"),nhbond,angle,dihedral)
        residue.coordinates = np.vstack((residue.coordinates,newh))
        residue.names = np.hstack((residue.names,["HZ3"]))
        residue.elements = np.hstack((residue.elements,["H"]))
        residue.natoms += 1

        if (ptm == "methylation") and (residue.name == "LYS"):
            newh = nerf(residue.find_coord("CD"),residue.find_coord("CE"),
                        residue.find_coord("NZ"),nhbond,109,150)
            residue.coordinates = np.vstack((residue.coordinates,newh))
            residue.names = np.hstack((residue.names,["HZ2"]))
            residue.elements = np.hstack((residue.elements,["H"]))
            residue.natoms += 1


    elif (residue.name in ["ARG"]) and (ptm != "dimethylation"):
        newh = nerf(residue.find_coord("NE"),residue.find_coord("CZ"),
                    residue.find_coord("NH1"),nhbond,120,0)
        residue.coordinates = np.vstack((residue.coordinates,newh))
        residue.names = np.hstack((residue.names,["H12"]))
        residue.elements = np.hstack((residue.elements,["H"]))
        residue.natoms += 1

        if ptm == "symmetric dimethylation":
            newh = nerf(residue.find_coord("NE"),residue.find_coord("CZ"),
                    residue.find_coord("NH2"),nhbond,120,0)
            residue.coordinates = np.vstack((residue.coordinates,newh))
            residue.names = np.hstack((residue.names,["H22"]))
            residue.elements = np.hstack((residue.elements,["H"]))
            residue.natoms += 1
        elif ptm == "asymmetric dimethylation":
            newh = nerf(residue.find_coord("NE"),residue.find_coord("CZ"),
                    residue.find_coord("NH2"),nhbond,120,180)
            residue.coordinates = np.vstack((residue.coordinates,newh))
            residue.names = np.hstack((residue.names,["H22"]))
            residue.elements = np.hstack((residue.elements,["H"]))
            residue.natoms += 1

    return

