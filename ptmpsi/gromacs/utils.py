def amber_to_gromacs_names(protein):
    for chain in protein.chains:
        for residue in chain.residues:
            names = []
            for index,name in enumerate(residue.names):
                string = "".join(name.split())
                if string == "HA2": names.append(" HA1"); continue
                if string == "HA3": names.append(" HA2"); continue
                if residue.name not in ["ALA"]:
                    if string == "HB2": names.append(" HB1"); continue
                    if string == "HB3": names.append(" HB2"); continue
                if string == "HG2": names.append(" HG1"); continue
                if string == "HG3": names.append(" HG2"); continue
                if residue.name not in ["PHE", "HIE", "HIP", "HID", "TYR", "ASH"]:
                    if string == "HD2": names.append(" HD1"); continue
                    if string == "HD3": names.append(" HD2"); continue
                if residue.name not in ["MET", "PHE", "TRP", "HIE", "TYR", "GLH" ]:
                    if string == "HE2": names.append(" HE1"); continue
                    if string == "HE3": names.append(" HE2"); continue
                if residue.name not in ["LYS", "TRP"]:
                    if string == "HZ2": names.append(" HZ1"); continue
                    if string == "HZ3": names.append(" HZ2"); continue
                if residue.name not in ["LEU","ASN"]: 
                    if string == "HD22": names.append("HD21"); continue
                    if string == "HD23": names.append("HD22"); continue
                if residue.name not in ["VAL"]: 
                    if string == "HG12": names.append("HG11"); continue
                    if string == "HG13": names.append("HG12"); continue
                if residue.name not in ["LEU"]:
                    if string == "HD11": names.append(" HD1"); continue
                    if string == "HD12": names.append(" HD2"); continue
                    if string == "HD13": names.append(" HD3"); continue
                if residue.name in ["ILE"]:
                    if string == "CD1": names.append(" CD"); continue
                names.append(residue.names[index])
            import numpy as np
            residue.names = np.array(names)

