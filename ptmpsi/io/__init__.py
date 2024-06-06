import numpy as np
from copy import deepcopy as copy
from ptmpsi.math import resdist

def writexyz(protein,xyzfile):
    with open(xyzfile,'w') as fh:
        fh.write(f"{protein.natoms}\n\n")
        for chain in protein.chains:
            for residue in chain.residues:
                for i in range(len(elements)):
                    fh.write("{}  {}  {}  {}\n".format(elements[i],*coordinates[i]))
    return
             

def writepdb(protein,pdbfile) :
    with open(pdbfile,'w') as fh:
        if protein.nssbonds > 0:
            for i,ssbond in enumerate(protein.ssbonds):
                fh.write("SSBOND {0: 3d} CYX {2} {1: 4d}    CYX {4} {3: 4d}{5: >43.2f}\n".format(i+1,*ssbond))
        i = 0
        for chain in protein.chains:
            for residue in chain.residues:
                for iatom in range(len(residue.names)):
                    i += 1
                    fh.write("ATOM  {:5d} {: >4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{: >24s}\n".format(i,residue.names[iatom],residue.name,residue.chain,residue.resid,*residue.coordinates[iatom],residue.elements[iatom]))
            fh.write("TER\n")
        fh.write("END")
    return

def writexyz(protein,xyzfile):
    with open(xyzfile,'w') as fh:
        print(protein.natoms, file=fh)
        print(" ", file=fh)
        for chain in protein.chains:
            for residue in chain.residues:
                for iatom in range(len(residue.elements)):
                    fh.write("{}  {: 16.8f}  {: 16.8f}  {: 16.8f}\n".format(residue.elements[iatom],*residue.coordinates[iatom]))
    return

def digestpdb(protein, interactive=False, delwat=True, delhet=True):
    protein.nonstandard = False
    protein.missing     = False
    nssbonds    = 0
    ssbonds     = []

    # Look for metadata
    for count,line in enumerate(protein.pdbfile):
        if line.split()[0] == "MODRES":
            protein.nonstandard = True
            print("\t A modified residue is present")
        if line.split()[0] == "SSBOND":
            nssbonds += 1
            ssbonds.append([int(line[17:21]),line[15:16],int(line[31:35]),line[29:30],float(line[73:78])])
        if line.split()[0] in ["ATOM","HETATM"]:
            start = count
            break

    resid = 0
    nmissing = 0
    natoms = 0
    alternate = False
    _chains = []
    _residues = []
    protein.nssbonds = nssbonds
    protein.ssbonds = copy(ssbonds)
    for line in protein.pdbfile[start:]:

        # Skip ANISOU lines
        if line[:6] == "ANISOU": continue

        # Skip CONECT lines
        if line[:6] == "CONECT": continue

        # Skip MASTER line
        if line[:6] == "MASTER": continue

        # Skip HETATM
        if (delhet) and (line[:6] == 'HETATM'): continue

        # Save last residue information
        if line.split()[0] in ["TER","END","ENDMDL"]:
            if resid > 0:
                residue = protein.savessbond(ssbonds,residue,resid,nmissing,chain)
                protein.addres(_residues,residue,atomid+1,names,elements,coordinates,chain,backbone)
                _residues, natoms, resid, nmissing = protein.addchain(_chains,chain,_residues,natoms,resid,nmissing)
            if line.split()[0] == "END": break
            continue
            
        if line[:4] == "ATOM" or line[:6] == "HETATM":
            # Check if we have a new chain
            if len(_residues) == 0 and resid == 0:
                chain   = line[21:22]
            # new chain was specified
            elif line[21:22] != chain:
                residue = protein.savessbond(ssbonds,residue,resid,nmissing,chain)
                protein.addres(_residues,residue,atomid+1,names,elements,coordinates,chain,backbone)
                _residues, natoms, resid, nmissing = protein.addchain(_chains,chain,_residues,natoms,resid,nmissing)
                chain   = line[21:22]
                    

            # Check if we have a new residue
            resnum  = int(line[22:26])
            if resnum != resid + nmissing:
                # Check if there are any missing internal residues
                tmpid = nmissing + resid + 1
                if tmpid != resnum:
                    protein.missing = True
                    if resid > 0:
                        print("\t Missing internal residues")
                    while tmpid != resnum:
                        tmpid += 1
                        nmissing += 1
                # Add information from previous residue
                if resid > 0:
                    residue = protein.savessbond(ssbonds,residue,resid,nmissing,chain)
                    protein.addres(_residues,residue,atomid+1,names,elements,coordinates,chain,backbone)
                # Increase residue counter and restart fields
                residue = line[17:20]
                resid += 1
                atomid = -1
                names = []
                elements = []
                coordinates = []
                alternate = False
                backbone = np.empty(3,dtype=int)


            # Read atom information
            atom    = line[12:16].strip()
            element = line[76:78].strip()

            # Check for alternate locations
            altcode = line[16:17]
            if alternate:
                if altcode != AorB:
                    alternate = False
                    continue
            elif altcode in ["A","B"]:
                alternate = True
                print("\t Warning: Atom {} of residue {} has alternate locations".format(atom,residue+str(resid+nmissing)))
                if interactive:
                    AorB = input("\t\t Select 'A' or 'B': ")
                    AorB = AorB.upper()
                    if AorB not in ["A","B"]:
                        raise ValueError("Did not understand alternate location: {}".format(AorB))
                else:
                    occupancy = float(line[54:60])
                    AorB = "A" if occupancy >= 0.5 else "B"
                    print("\t\t Selecting location '{}' with occupancy {:6.2f}".format(AorB,max(occupancy,1-occupancy)))
                if altcode != AorB:
                    continue

            # Increase counters
            atomid += 1
            natoms += 1

            # Append elements
            names.append(atom)
            elements.append(element)
            coordinates.append([float(line[30:38]), 
                                float(line[38:46]), 
                                float(line[46:54])])

            # Define backbone
            if atom == 'N':
                backbone[0] = atomid
            elif atom == 'CA':
                backbone[1] = atomid
            elif atom == 'C':
                backbone[2] = atomid


    # Save last residue information in case PDB file ended abruptly
    if len(_residues) > 0:
        residue = protein.savessbond(ssbonds,residue,resid,nmissing,chain)
        protein.addres(_residues,residue,atomid+1,names,elements,coordinates,chain,backbone)
        _residues,natoms,resid,nmissing = protein.addchain(_chains,chain,_residues,natoms,resid,nmissing)
    protein.chains = copy(_chains)
    protein.nchains = len(protein.chains)

    # Remove water residues
    if delwat: protein.delwaters()

    # Update number of residues and atoms
    protein.update()

    # Check SSBOND distance
    for issbond in protein.ssbonds:
        for chain in protein.chains:
            if chain.name == issbond[1]:
                residue1 = chain.residues[issbond[0]-1]
            if chain.name == issbond[3]:
                residue2 = chain.residues[issbond[2]-1]
        _distance,atom1,atom2 = resdist(residue1,residue2)
        if np.isclose(_distance,issbond[4],0.1): continue
        print("\t Warning!!!")
        print("\t\t Expected SSBOND with {:4.2f} A bond length, but got {:4.2f} A instead".format(issbond[4],_distance))


    return
