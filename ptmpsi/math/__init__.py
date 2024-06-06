import numpy as np
import math
from copy import deepcopy as copy
from ptmpsi.constants import amidebond, nhbond, amideangle

aminolist = [
            "ACE",
            "ALA",
            "ARG",
            "ASH",
            "ASN",
            "ASP",
            "CYM",
            "CYS",
            "CYX",
            "GLH",
            "GLN",
            "GLU",
            "GLY",
            "HID",
            "HIE",
            "HIP",
            "HIS",
            "ILE",
            "LEU",
            "LYN",
            "LYS",
            "MET",
            "NHE",
            "NME",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL"
        ]
####

def norm(vec):
    return math.sqrt( vec[0]**2 + vec[1]**2 + vec[2]**2 )

def resdist(residue1,residue2):
    """
    Computes the minimum distance between two residues
    """
    distance = np.finfo(float).max
    for iatom,icoor in enumerate(residue1.coordinates):
        for jatom,jcoor in enumerate(residue2.coordinates):
            _distance = norm(icoor-jcoor)
            if _distance < distance:
                distance = _distance
                ipos = iatom
                jpos = jatom
    return distance,ipos,jpos

####

def nerf(atoma,atomb,atomc,bond,angle,torsion):
    """
    Obtain cartesian coordinates for a fourth atom
    based on the bond length CD, angle BCD, and
    torsion ABCD, using the Natural Extension Reference Frame
    method.
    """
    _angle = np.pi - np.radians(angle)
    _torsion = np.radians(torsion)
    sinang = np.sin(_angle)
    atomd = bond*np.array([ np.cos(_angle), 
        sinang*np.cos(_torsion), 
        sinang*np.sin(_torsion) ])

    AB = atomb - atoma
    BC = atomc - atomb
    bc = BC/norm(BC)
    n  = np.cross(AB,bc)
    n  = n/norm(n)
    nxbc = np.cross(n,bc)

    R = np.column_stack((bc,nxbc,n))

    return np.dot(R,atomd) + atomc




def appendc(chain,residue,psi=None):
    """
    Attach residue to C-terminus of chain
    """

    if len(chain.residues) == 0:
        coords = residue.coordinates
        return coords

    # Find C-terminus
    cterminus = chain.residues[-1]
    if cterminus.name in ["NME","NHE"]:
        print("\t Warning: a new residue cannot be appended to '{}'".format(cterminus.name))
        print("\t\t Nothing will be done!")
        return None

    # Find Alpha-carbon and carbonyl group
    ca = cterminus.find_coord("CA")
    c  = cterminus.find_coord("C")
    o  = cterminus.find_coord("O")

    if cterminus.name in ["ABA", "EAC"]:
        n = cterminus.find_coord("CB")
    else:
        n = cterminus.find_coord("N")

    # If C-terminus had OXT oxygen, take its position to attach the next residue
    # the C atom will be placed at the origin
    if np.any(cterminus.names == "OXT"):
        vector1 = cterminus.find_coord("OXT") - c
        vector1 = vector1/norm(vector1)*amidebond
        mask = [*range(0,len(cterminus.names))]
        mask.pop(cterminus.find("OXT"))
        cterminus.coordinates = cterminus.coordinates[mask]
        cterminus.names = cterminus.names[mask]
        cterminus.elements = cterminus.elements[mask]
        cterminus.natoms = len(cterminus.names)

        ca = cterminus.find_coord("CA")
        c  = cterminus.find_coord("C")
        o  = cterminus.find_coord("O")
        if cterminus.name in ["ABA", "EAC"]:
            n = cterminus.find_coord("CB")
        else:
            n = cterminus.find_coord("N")

    # Otherwise, bisect the angle CA-C-O to determine position
    else:
        _c  = c - c
        _ca = ca - c
        _o  = o - c
        bisector = norm(_o)*_ca + norm(_ca)*_o
        vector1 = -amidebond*bisector/norm(bisector)

    # Attachment position should be already in the template
    # The N atom will be placed at the origin
    vector2 = residue.nattach - residue.find_coord("N")

    # Obtain rotation matrix and ensure antiparallel arrangement
    R = rotmatvec(vector1,-vector2)

    # Rotate coordinates and translate to original frame
    newcoords = np.dot(residue.coordinates-residue.find_coord("N"),R.T) + c + vector1

    # Check placement of H atom. If needed, rotate along the amide bond to generate
    # anti H-N-C-O geometry
    if residue.name not in ["PRO"]:
        dihedral = get_torsion(newcoords[residue.find("H")],
                newcoords[residue.find("N")],c,o)
        if abs(dihedral) < 175.0:
            shift = copy(newcoords[residue.find("N")])
            k = -(c - shift)/norm(c-shift)
            dihedral = np.radians(180.0-dihedral)
            for i,coord in enumerate(newcoords):
                if i == residue.find("N"): continue
                coord = np.cos(dihedral)*(coord - shift) + (1-np.cos(dihedral))*k*np.dot(coord-shift,k) + np.cross(k,coord-shift)*np.sin(dihedral) + shift
                newcoords[i] = coord

    if psi is not None:

        # Check Psi angle
        dihedral = get_torsion(n, ca, c,
            newcoords[residue.find("N")]) - psi

        # Rotate, if necessary
        if np.abs(dihedral) > 0.0:
            dihedral = -np.radians(dihedral)
            R = rotmataxis(cterminus.find_coord("C")-cterminus.find_coord("CA"),dihedral)
            newcoords = np.dot(newcoords-cterminus.find_coord("CA"),R.T) + cterminus.find_coord("CA")
            cterminus.coordinates[cterminus.find("O")] = np.dot(R,cterminus.find_coord("O")-cterminus.find_coord("CA")) + cterminus.find_coord("CA")


    return newcoords


def get_torsion(atoma,atomb,atomc,atomd):
    b0 = atoma - atomb
    b1 = atomc - atomb
    b2 = atomd - atomc
    b1 /= norm(b1)
    v = b0 - np.dot(b0,b1)*b1
    w = b2 - np.dot(b2,b1)*b1
    x = np.dot(v,w)
    y = np.dot(np.cross(b1,v),w)
    return np.degrees(np.arctan2(y,x))


def prependn(chain,residue,phi=None):
    """
    Attach residue to N-terminus
    """

    if len(chain.residues) == 0:
        return residue.coordinates

    # Find N-terminus
    nterminus = chain.residues[0]
    if nterminus.name in ["ACE"]:
        print("\t Warning: A new residue cannot be prepended to '{}'".format(nterminus.name))
        print("\t\t Nothing will be done!")
        return None

    # Find backbone
    if nterminus.name in ["EAC"]:
        n  = copy(nterminus.find_coord("N"))
        ca = copy(nterminus.find_coord("CE"))
        c  = copy(nterminus.find_coord("CD"))
    elif nterminus.name in ["ABA"]:
        n  = copy(nterminus.find_coord("N"))
        ca = copy(nterminus.find_coord("CG"))
        c  = copy(nterminus.find_coord("CB"))
    else:
        n  = copy(nterminus.find_coord("N"))
        ca = copy(nterminus.find_coord("CA"))
        c  = copy(nterminus.find_coord("C"))

    # Check if N-terminus was protonated
    mask = []
    for i,atom in enumerate(nterminus.names):
        if atom in ["H1","H2","H3"]: continue
        mask.append(i)

    # Remove hydrogens from lists
    if len(mask) < len(nterminus.names):
        nterminus.names = nterminus.names[mask]
        nterminus.coordinates = nterminus.coordinates[mask]
        nterminus.elements = nterminus.elements[mask]
        nterminus.natoms = len(nterminus.names)
    try:
        hpos = nterminus.find("H")
    except:
        hpos = nterminus.natoms
        nterminus.natoms += 1

        _temp = np.empty(nterminus.natoms, dtype='U4')
        _temp[:nterminus.natoms-1] = nterminus.names.copy()
        _temp[-1] = "H"
        nterminus.names = _temp.copy()

        _temp = np.empty(nterminus.natoms, dtype='U4')
        _temp[:nterminus.natoms-1] = nterminus.elements.copy()
        _temp[-1] = "H"
        nterminus.elements = _temp.copy()
        
        _temp = np.empty((nterminus.natoms,3), dtype=float)
        _temp[:nterminus.natoms-1] = nterminus.coordinates.copy()
        _temp[-1] = np.array([0.0, 0.0 ,0.0])
        nterminus.coordinates = _temp.copy()

    # Determine attachment vector for N-terminus
    vector1 = nerf(c,ca,n,amidebond,amideangle,phi) - n

    # Determine attachment vector for new residue
    vector2 = residue.cattach - residue.find_coord("C")

    # Obtain rotation matrix and ensure antiparallel arrangement
    R = rotmatvec(vector1,-vector2)

    # Rotate coordinates and translate to original frame
    newcoords = np.dot(residue.coordinates-residue.find_coord("C"),R.T) + n + vector1

    # Ensure that the amide bond is planar
    dihedral = get_torsion(newcoords[residue.find("O")],
            newcoords[residue.find("C")],n,ca)

    if abs(dihedral) > 5.0:
        dihedral = np.radians(dihedral)
        R = rotmataxis(n-newcoords[residue.find("C")],dihedral)
        temp = newcoords - newcoords[residue.find("C")]
        newcoords = np.dot(temp,R.T) + newcoords[residue.find("C")]

    # Add hydrogen antiparallel to carbonyl
    vector1 = nerf(newcoords[residue.find("O")],newcoords[residue.find("C")],n,nhbond,120,180)
    nterminus.coordinates[hpos] = copy(vector1)
    nterminus.backbone[0] = nterminus.find("N")
    nterminus.backbone[1] = nterminus.find("CA")
    nterminus.backbone[2] = nterminus.find("C")


    if phi is not None:
        # Check Phi angle
        dihedral = get_torsion(newcoords[residue.find("C")],
            n, ca, c) - phi

        # Rotate, if necessary
        if np.abs(dihedral) > 0.5:
            dihedral = -np.radians(dihedral)
            R = rotmataxis(nterminus.find_coord("N")-nterminus.find_coord("CA"),dihedral)
            newcoords = np.dot(newcoords-nterminus.find_coord("CA"),R.T) + nterminus.find_coord("CA")
            nterminus.coordinates[nterminus.find("H")] = np.dot(R,nterminus.find_coord("H")-nterminus.find_coord("CA")) + nterminus.find_coord("CA")

        dihedral = get_torsion(newcoords[residue.find("C")],
            nterminus.find_coord("N"),
            nterminus.find_coord("CA"),
            nterminus.find_coord("C"))


    return newcoords


def rotmatvec(vector1,vector2):
    """
    Return the rotation matrix to align vector2
    into vector1
    """
    vhat1 = vector1/norm(vector1)
    vhat2 = vector2/norm(vector2)
    cross = np.cross(vhat2,vhat1)
    sin = norm(cross)
    cos = np.dot(vhat2,vhat1)
    if np.isclose(cos,1,1.0E-6):
        R = np.eye(3)
    elif np.isclose(cos,-1,1.0E-6):
        R = -np.eye(3)
    else:
        skew = np.zeros((3,3))
        skew[0] = np.array([0,-cross[2],cross[1]])
        skew[1] = np.array([cross[2],0,-cross[0]])
        skew[2] = np.array([-cross[1],cross[0],0])
        R = np.eye(3) + skew + np.dot(skew,skew)/(1+cos)
    return R


def rotmataxis(k,theta):
    """
    Obtain the rotation matrix through an angle theta 
    counterclockwise about the axis k
    """
    khat = k/norm(k)
    K = np.zeros((3,3))
    K[0] = [0,-khat[2],khat[1]]
    K[1] = [khat[2],0,-khat[0]]
    K[2] = [-khat[1],khat[0],0]
    return np.eye(3) + np.sin(theta)*K + (1-np.cos(theta))*K@K


def alignres(residue1,residue2):
    """
    Rotates and translates residue 2 to align as 
    best as possible the backbone atoms CA, C, and N
    with residue 1.
    """

    # Get backbone coordinates
    refcoords = residue1.coordinates[ residue1.backbone ]
    refcoords = np.vstack([refcoords, residue1.find_coord("O")])
    tocoords  = residue2.coordinates[ residue2.backbone ]
    tocoords = np.vstack([tocoords, residue2.find_coord("O")])

    # Set CA at the origin
    refshift = copy(refcoords[1])
    refcoords -= refshift

    toshift = copy(tocoords[1])
    tocoords  -= toshift

    # Get PHI torsion angles
    refphi = get_torsion(refcoords[0],refcoords[1],refcoords[2],refcoords[3])
    tophi = get_torsion(tocoords[0],tocoords[1],tocoords[2],tocoords[3])

    # Rotate over PHI torsion to match reference
    angle = -np.radians(refphi - tophi)
    if abs(angle) > 0.1:
        Rphi = rotmataxis(tocoords[2],angle)
        tocoords[0] = np.dot(Rphi,tocoords[0])


    # Get Covariance Matrix
    covar = np.dot(tocoords[:3].T, refcoords[:3])

    # Build F matrix
    F = np.array([
        [covar[0,0]+covar[1,1]+covar[2,2],covar[1,2]-covar[2,1],covar[2,0]-covar[0,2],covar[0,1]-covar[1,0]],
        [covar[1,2]-covar[2,1],covar[0,0]-covar[1,1]-covar[2,2],covar[0,1]+covar[1,0],covar[2,0]+covar[0,2]],
        [covar[2,0]-covar[0,2],covar[0,1]+covar[1,0],-covar[0,0]+covar[1,1]-covar[2,2],covar[1,2]+covar[2,1]],
        [covar[0,1]-covar[1,0],covar[0,2]+covar[2,0],covar[1,2]+covar[2,1],-covar[0,0]-covar[1,1]+covar[2,2]]
    ])

    # Eigendecomposition of F matrix
    l, qs = np.linalg.eigh(F)

    # Select Eigenvector corresponding to the largest Eigenvalue
    q = qs[:,-1]

    # Build rotation matrix
    R = np.array([
        [q[0]**2+q[1]**2-q[2]**2-q[3]**2,2*(q[1]*q[2]-q[0]*q[3]),2*(q[1]*q[3]+q[0]*q[2])],
        [2*(q[1]*q[2]+q[0]*q[3]),q[0]**2-q[1]**2+q[2]**2-q[3]**2,2*(q[2]*q[3]-q[0]*q[1])],
        [2*(q[1]*q[3]-q[0]*q[2]),2*(q[2]*q[3]+q[0]*q[1]),q[0]**2-q[1]**2-q[2]**2+q[3]**2]
        ])

    # Rotate coordinates
    newcoords = residue2.coordinates - toshift
    if abs(angle) > 0.1:
        for i,coord in enumerate(newcoords):
            if residue2.elements[i,1] in ["CA","C","O"]: continue
            newcoords[i] = np.dot(Rphi,coord)
    newcoords = np.dot( newcoords, R.T)


    # Translate coordinates to reference ones
    newcoords += refshift

    return newcoords



def rotate_chi1(residue,atoms,chi1):
    if (atoms is None) or (chi1 == 0): return
    tmp = residue.coordinates[atoms] - residue.coordinates[atoms[1]]
    angle = -np.radians(chi1)
    R = rotmataxis(tmp[2],angle)
    mask = []
    for i,name in enumerate(residue.names):
        if name in ["N","H","CA","HA","C","O"]: continue
        mask.append(i)
    residue.coordinates[mask] = np.dot(residue.coordinates[mask]-residue.coordinates[atoms[1]],R.T) + residue.coordinates[atoms[1]]

    internal = False
    for iatom in range(len(residue.coordinates)):
        for jatom in range(iatom+1,len(residue.coordinates)):
            if norm(residue.coordinates[iatom] - residue.coordinates[jatom]) < 1.0:
                internal = True
                break
    return internal


def rotate_chi2(residue,atoms,chi2):
    if (atoms is None) or (chi2 == 0): return
    tmp = residue.coordinates[atoms] - residue.coordinates[atoms[1]]
    angle = -np.radians(chi2)
    R = rotmataxis(tmp[2],angle)
    mask = []
    for i,name in enumerate(residue.names):
        if name in ["N","H","CA","HA","C","O","CB","HB","HB1","HB2","HB3"]: continue
        if (residue.name == "ILE") and (name in ["CG2","HG21","HG22","HG23"]): continue
        mask.append(i)
    residue.coordinates[mask] = np.dot(residue.coordinates[mask]-residue.coordinates[atoms[1]],R.T) + residue.coordinates[atoms[1]]

    internal = False
    for iatom in range(len(residue.coordinates)):
        for jatom in range(iatom+1,len(residue.coordinates)):
            if norm(residue.coordinates[iatom] - residue.coordinates[jatom]) < 1.0:
                internal = True
                break

    return internal


def find_clashes(proteins):
    _all = []
    nclashes = 0

    for protein in proteins:
        for chain in protein.chains:
            _all += chain.residues

    for iresidue in _all:
        nclashes += find_clashes_residue(iresidue,proteins,allres=_all)

    return



def find_clashes_residue(iresidue,proteins,printing=True,allres=None):

    _all = [] if allres is None else allres

    _proteins = []
    _iprotein = -1
    nclashes = 0
    for protein in proteins:
        _iprotein += 1
        for chain in protein.chains:
            if allres is None: _all += chain.residues
            for residue in chain.residues:
                _proteins.append( _iprotein )

    for jres in range(len(_all)):
        jresidue = _all[jres]
        if jresidue == iresidue: continue
        ires = jres

        # Test if residues are contiguous
        _samechain = iresidue.chain == jresidue.chain
        _contiguous = _samechain and (abs(iresidue.resid-jresidue.resid) == 1)

        # Random distance
        _distance = norm(iresidue.coordinates[0]-jresidue.coordinates[0])

        # If large, it is safe to assume no clashes between residues
        if _distance > 10.0: continue

        for iatom,iname in enumerate(iresidue.names):
            _sg = iname == "SG"
            _n  = iname in ["N","H"]
            _c  = iname in ["C","O"] 
            _ih = iname[0] == "H"
            for jatom,jname in enumerate(jresidue.names):
                _ssbond = _sg and jname == "SG"
                _amide  = _contiguous and ((_n and jname in ["C","O","CA"]) or (_c and jname in ["N","H","CA"]))
                _jh     = jname[0] == "H"

                # Skip amide bonds and possible ssbonds
                if _amide or _ssbond: continue

                # Get distance between atoms
                _distance = norm(iresidue.coordinates[iatom]-jresidue.coordinates[jatom])
                # Heavy atoms, longer threshold
                threshold = 1.4 if (_ih or _jh) else 2.3

                if _distance < threshold:
                    nclashes += 1
                    if not printing: continue
                    print("\t Warning: Possible clash, distance: {:8.3f} A".format(_distance))
                    print("\t\t Atom {}:{}:{}{}:{} and {}:{}:{}{}:{}".format(
                        _proteins[ires]+1,iresidue.chain,iresidue.name,str(iresidue.resid),iresidue.names[iatom],
                        _proteins[jres]+1,jresidue.chain,jresidue.name,str(jresidue.resid),jresidue.names[jatom]))

    return nclashes
