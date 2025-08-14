# selection_features.py

import os
import json

import numpy as np
import pandas as pd
import tempfile
from tabulate import tabulate
from joblib import Parallel, delayed
import MDAnalysis as mda

# bring in only the constants & MDAnalysis classes you need
from pca_config import Universe, Writer, distances, mda_rmsd, STANDARD_NAMES

import mdtraj as md


def get_unique(ag):
    return ag.unique() if callable(ag.unique) else ag.unique

def write_pretty_tsv(df: pd.DataFrame, fn: str):
    """
    Write a tab-separated file with padded, aligned columns.
    """
    tsv_text = tabulate(
        df.values,
        headers=df.columns,
        tablefmt="tsv",
        floatfmt=".6f",
        colalign=("left",) + ("right",)*(df.shape[1]-1),
    )
    with open(fn, "w") as f:
        f.write(tsv_text)
    print(f"TSV written to {fn}")

def write_pretty_csv(df: pd.DataFrame, fn: str):
    """
    Write a comma-separated file with padded, aligned columns.
    """
    # Generate a “TSV” block with padded widths, then replace tabs with commas.
    tsv_text = tabulate(
        df.values,
        headers=df.columns,
        tablefmt="tsv",
        floatfmt=".6f",
        colalign=("left",) + ("right",)*(df.shape[1]-1),
    )
    csv_text = tsv_text.replace("\t", ",")
    with open(fn, "w") as f:
        f.write(csv_text)
    print(f"CSV written to {fn}")

def load_and_preprocess_subsystem_mda(system_name, system_info, subsample_step):
    """
    Load a single MDAnalysis Universe and subsample every `subsample_step`‑th frame.

    Returns:
      - u: the aligned MDAnalysis Universe,
      - frame_indices: list of subsampled frame indices.
    """
    top_file = system_info["top"]
    traj_files = system_info["traj_list"]
    ptm_type = system_info.get("ptm_type", "unknown")

    print(f"\n=== Loading system '{system_name}' ===")
    print(f"    PTM Type: {ptm_type}, #runs={len(traj_files)}")

    # Load the Universe (topology + all trajectory files)
    u = mda.Universe(top_file, *traj_files)
    u.topology_filename = top_file
    n_frames = len(u.trajectory)
    print(f"    Total frames in Universe: {n_frames}")

    # Print segment IDs for information
    print("Segment IDs in the topology:")
    for seg in u.segments:
        print(f"    segid: {seg.segid}")

    # Align the trajectory to its first frame using protein atoms.
    # u.trajectory[0]  # set the reference frame
    # alignment = align.AlignTraj(u, u, select="name N CA C", in_memory=True)
    # alignment.run()

    # Subsample frame indices manually
    frame_indices = list(range(0, n_frames, subsample_step))
    print(f"    Subsampling with step={subsample_step} => {len(frame_indices)} frames used.")

    return u, frame_indices

def select_atoms_for_analysis(
    u,
    ptm_resnames=None,
    standard_resnames=STANDARD_NAMES,
    gap_segids=None,
    prk_segids=None,
    cp12_segids=None,
    nad_resnames=None,
    radius_ptm=3.0,
    radius_nad=3.0,
    radius_interfaces=2.0
):


    if ptm_resnames is None:
        ptm_resnames = "AUTO"
    if gap_segids is None:
        gap_segids = ["A","D","E","F","I","K","N","O"]
    if prk_segids is None:
        prk_segids = ["B","G","J","L"]
    if cp12_segids is None:
        cp12_segids = ["C","H","M","P"] 
    if nad_resnames is None:
        nad_resnames = ["NAD","NADH","NADPH","NAI"]

    # If ptm_resnames is "AUTO", detect non-standard residues.
    if isinstance(ptm_resnames, str) and ptm_resnames.upper() == "AUTO":
        auto_list = []
        for res in u.residues:
            if res.resname.upper() not in standard_resnames:
                auto_list.append(res.resname.upper())
        ptm_resnames = list(set(auto_list))
        print("Auto-detected non-standard PTM resnames:", ptm_resnames)
    
    ptm_group = u.atoms[[]]
    for rname in ptm_resnames:
        ptm_group += u.select_atoms(f"resname {rname}")
    ptm_neighbors = u.select_atoms(f"byres (around {radius_ptm} group ptm_group) and not resname HOH WAT TIP3 SOL MG ZN NA CL K CA", ptm_group=ptm_group)
    ptm_total = ptm_group + ptm_neighbors
    
    nad_group = u.atoms[[]]
    for rname in nad_resnames:
        nad_group += u.select_atoms(f"resname {rname}")
    nad_neighbors = u.select_atoms(f"byres (around {radius_nad} group nad_group) and not resname HOH WAT TIP3 SOL MG ZN NA CL K CA", nad_group=nad_group)
    nad_total = nad_group + nad_neighbors
    
    def segid_to_query(segids):
        return " or ".join([f"segid {s}" for s in segids])
    gap_sel  = u.select_atoms(segid_to_query(gap_segids))
    prk_sel  = u.select_atoms(segid_to_query(prk_segids))
    cp12_sel = u.select_atoms(segid_to_query(cp12_segids))
    
    gap_around_prk = u.select_atoms(f"byres (around {radius_interfaces} group gap_sel) and not resname HOH WAT TIP3 SOL MG ZN NA CL K CA", gap_sel=gap_sel) & prk_sel
    prk_around_gap = u.select_atoms(f"byres (around {radius_interfaces} group prk_sel) and not resname HOH WAT TIP3 SOL MG ZN NA CL K CA", prk_sel=prk_sel) & gap_sel
    gap_prk_if = gap_around_prk + prk_around_gap
    
    prk_around_cp12 = u.select_atoms(f"byres (around {radius_interfaces} group prk_sel) and not resname HOH WAT TIP3 SOL MG ZN NA CL K CA", prk_sel=prk_sel) & cp12_sel
    cp12_around_prk = u.select_atoms(f"byres (around {radius_interfaces} group cp12_sel) and not resname HOH WAT TIP3 SOL MG ZN NA CL K CA", cp12_sel=cp12_sel) & prk_sel
    prk_cp12_if = prk_around_cp12 + cp12_around_prk
    
    gap_around_cp12 = u.select_atoms(f"byres (around {radius_interfaces} group gap_sel) and not resname HOH WAT TIP3 SOL MG ZN NA CL K CA", gap_sel=gap_sel) & cp12_sel
    cp12_around_gap = u.select_atoms(f"byres (around {radius_interfaces} group cp12_sel) and not resname HOH WAT TIP3 SOL MG ZN NA CL K CA", cp12_sel=cp12_sel) & gap_sel
    gap_cp12_if = gap_around_cp12 + cp12_around_gap

    interfaces = {
        "gap_prk": get_unique(gap_prk_if),
        "prk_cp12": get_unique(prk_cp12_if),
        "gap_cp12": get_unique(gap_cp12_if)
    }
    
    final_selection = ptm_total + nad_total + prk_cp12_if + gap_cp12_if
    try:
        final_unique = final_selection.unique()
    except TypeError:
        final_unique = final_selection.unique
    print(f"  PTM raw atoms          : {ptm_group.n_atoms}")
    print(f"  PTM neighbors (3Å)    : {ptm_neighbors.n_atoms}")
    print(f"  PTM total (group+shell): {(ptm_group + ptm_neighbors).n_atoms}")
    return {"final": final_unique, "interfaces": interfaces, "ptm_total": ptm_group + ptm_neighbors, "nad_total": nad_group + nad_neighbors}

def extract_and_save_features(
    system_name,
    info,
    step,
    feat_dir,
    radius_ptm,
    radius_nad,
    radius_interfaces,
    dssp_window,
    contact_cutoff
):
    u, frames = load_and_preprocess_subsystem_mda(system_name, info, step)
    sel = select_atoms_for_analysis(
        u,
        ptm_resnames      = info.get("ptm_resnames","AUTO"),
        radius_ptm        = radius_ptm,
        radius_nad        = radius_nad,
        radius_interfaces = radius_interfaces
    )
    if sel["final"].n_atoms == 0:
        # no atoms → empty features _and_ empty PTM outputs
#        feats      = np.zeros((len(frames), 0))
#        names      = []
#        ptm_mat    = np.zeros((len(frames), 0))
#        ptm_labels = []
        feats      = np.zeros((len(frames), 1))
        names      = ["ZeroFeature"]
        ptm_mat    = np.zeros((len(frames), 1))
        ptm_labels = ["ZeroFeature"]
    else:
        feats, names, ptm_mat, ptm_labels = extract_features_mda(
            u,
            sel["final"],
            frames,
            compute_pairwise_dist=False,
            native_contact_cutoff=contact_cutoff,
            dssp_window=dssp_window,
            ptm_resnames_for_contacts=info.get("ptm_resnames",[])
        )
    os.makedirs(feat_dir, exist_ok=True)
    np.save(os.path.join(feat_dir, f"{system_name}_feats.npy"), feats)
    with open(os.path.join(feat_dir, f"{system_name}_featnames.json"), "w") as f:
        json.dump(names, f)
    # NEW: save PTM contacts & labels
    np.save(os.path.join(feat_dir, f"{system_name}_ptm_matrix.npy"), ptm_mat)
    with open(os.path.join(feat_dir, f"{system_name}_ptm_labels.json"), "w") as f:
        json.dump(ptm_labels, f)
    return system_name

def compute_rg(coords):
    com = coords.mean(axis=0)
    sq = ((coords - com)**2).sum(axis=1)
    return np.sqrt(sq.mean()) * 0.1

def calc_dihedral(p0, p1, p2, p3):
    """
    Compute dihedral angle (in radians) defined by four positions.
    """
    b0 = p1 - p0
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.arctan2(y, x)

def compute_avg_phi_psi(u, residues):
    """
    Compute average φ and ψ angles (in radians) over a sorted list of protein residues.
    For residue i: φ = dihedral(C(i-1), N(i), CA(i), C(i)) and
    ψ = dihedral(N(i), CA(i), C(i), N(i+1)).
    Returns (avg_phi, avg_psi).
    """
    phis = []
    psis = []
    res_list = sorted(residues, key=lambda r: (r.segid, r.resid))
    for i in range(1, len(res_list)-1):
        try:
            C_prev = res_list[i-1].atoms.select_atoms("name C")[0].position
            N_i    = res_list[i].atoms.select_atoms("name N")[0].position
            CA_i   = res_list[i].atoms.select_atoms("name CA")[0].position
            C_i    = res_list[i].atoms.select_atoms("name C")[0].position
            N_next = res_list[i+1].atoms.select_atoms("name N")[0].position
            phi = calc_dihedral(C_prev, N_i, CA_i, C_i)
            psi = calc_dihedral(N_i, CA_i, C_i, N_next)
            phis.append(phi)
            psis.append(psi)
        except Exception as e:
            print(f"Error computing dihedral for residue index {i}: {e}")
            continue
    if len(phis) == 0:
        return 0.0, 0.0
    return np.mean(phis), np.mean(psis)

def extract_features_mda(
    u,
    sel_atoms,
    frame_indices,
    compute_pairwise_dist=False,
    native_contact_cutoff=3.0,
    dssp_window=25,
    ptm_resnames_for_contacts=None
):
    """
    Extract features from selected atoms over specified frames.

    Features include:
      - Optional pairwise distances (flattened upper triangle).
      - Radius of gyration (Rg).
      - RMSD relative to the first frame.
      - Sliding-window B-factor (frame-dependent).
      - Manual average φ and ψ (computed from sel_atoms residues).
      - DSSP secondary structure fractions (via MDTraj).
      - Native contacts per PTM site.
    """
    n_frames = len(frame_indices)
    n_sel    = sel_atoms.n_atoms
    print(f"Extracting features from {n_sel} selected atoms, frames={n_frames}.")

    # Optional pairwise distances
    pair_indices = []
    if compute_pairwise_dist:
        for i in range(n_sel):
            for j in range(i+1, n_sel):
                pair_indices.append((i, j))
    n_pairs = len(pair_indices)
    dist_arr = np.zeros((n_frames, n_pairs)) if compute_pairwise_dist and n_pairs>0 else None

    # Prepare arrays
    rg_arr   = np.zeros(n_frames)
    rmsd_arr = np.zeros(n_frames)

    # Reference coordinates for RMSD
    u.trajectory[frame_indices[0]]
    ref_coords = sel_atoms.positions.copy()

    # ---- COM per residue for B-factor window ----
    coms_dict = {}
    for res in sel_atoms.residues:
        label = f"{res.resname}-{res.segid}-{res.resid}"
        coords = []
        for fi in frame_indices:
            u.trajectory[fi]
            coords.append(res.atoms.center_of_mass())
        coms_dict[label] = np.array(coords)  # shape (n_frames, 3)

    # ---- Sliding-window B-factor per frame ----
    window_half = dssp_window
    b_factors = np.zeros(n_frames)
    for i in range(n_frames):
        start = max(0, i - window_half)
        end   = min(n_frames, i + window_half + 1)
        window_means = []
        for coords in coms_dict.values():
            segment = coords[start:end]
            mean_pos = segment.mean(axis=0)
            rmsf = np.sqrt(np.mean(np.sum((segment - mean_pos)**2, axis=1)))
            
            # B-factor in Å² → convert to nm² by multiplying by (0.1)² = 0.01
            bval_A2  = (rmsf**2) * (8.0 * np.pi**2) / 3.0
            bval_nm2 = bval_A2 * 0.01
            window_means.append(bval_nm2)
        b_factors[i] = np.mean(window_means)

    # ---- Manual average φ/ψ ----
    prot_residues = [res for res in sel_atoms.residues if res.atoms.select_atoms("name N CA C").n_atoms == 3]
    avg_phi = np.zeros(n_frames)
    avg_psi = np.zeros(n_frames)
    for i, fi in enumerate(frame_indices):
        u.trajectory[fi]
        phi, psi = compute_avg_phi_psi(u, prot_residues)
        avg_phi[i] = phi
        avg_psi[i] = psi

    # ---- DSSP via MDTraj (one frame at a time) ----
    frac_helix = np.zeros(n_frames)
    frac_sheet = np.zeros(n_frames)
    frac_coil  = np.zeros(n_frames)
    # prepare DSSP topology once
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp:
        topo_pdb = tmp.name
        with mda.Writer(topo_pdb, sel_atoms.n_atoms) as W:
            u.trajectory[frame_indices[0]]
            W.write(sel_atoms)
    md_top = md.load_pdb(topo_pdb).topology
    os.remove(topo_pdb)
    helix_set = {"H","G","I"}
    sheet_set = {"E","B"}
    # compute DSSP per frame
    for idx, fi in enumerate(frame_indices):
        u.trajectory[fi]
        coords = sel_atoms.positions / 10.0  # nm
        one = md.Trajectory(coords[np.newaxis, :, :], md_top)
        ss = md.compute_dssp(one)[0]
        frac_helix[idx] = np.mean([c in helix_set for c in ss])
        frac_sheet[idx] = np.mean([c in sheet_set for c in ss])
        frac_coil[idx]  = 1.0 - frac_helix[idx] - frac_sheet[idx]

    # ---- Native contacts per PTM site ----
    auto_flag = isinstance(ptm_resnames_for_contacts, str) and ptm_resnames_for_contacts.upper()=="AUTO"
    ptm_sites = []
    if auto_flag:
        for res in u.residues:
            if res.resname.upper() not in STANDARD_NAMES:
                ptm_sites.append(res)
    else:
        for res in u.residues:
            if res.resname.upper() in (ptm_resnames_for_contacts or []):
                ptm_sites.append(res)
    n_sites = len(ptm_sites)
    ptm_contact_matrix = np.zeros((n_frames, n_sites), dtype=float)
    ptm_site_labels   = [f"{r.resname}-{r.segid}-{r.resid}" for r in ptm_sites]
    for i_site, res in enumerate(ptm_sites):
        for idx, fi in enumerate(frame_indices):
            u.trajectory[fi]
            ptm_atoms = res.atoms
            prot = u.select_atoms("protein") - ptm_atoms
            contacts = set()
            for other in prot.residues:
                d = distances.distance_array(ptm_atoms.positions,
                                             other.atoms.positions,
                                             box=u.trajectory.ts.dimensions)
                if np.any(d < native_contact_cutoff):
                    contacts.add((other.segid, other.resid))
            ptm_contact_matrix[idx, i_site] = len(contacts)

    # ---- Rg, RMSD, Pairwise distances ----
    for idx, fi in enumerate(frame_indices):
        u.trajectory[fi]
        pos = sel_atoms.positions.copy()
        rg_arr[idx]   = compute_rg(pos)
        # now in nm        
        rmsd_arr[idx] = mda_rmsd(pos, ref_coords, center=True, superposition=True) * 0.1
        if dist_arr is not None:
            dmat = distances.distance_array(pos, pos, box=u.trajectory.ts.dimensions)
            dist_arr[idx, :] = [dmat[i, j] for (i, j) in pair_indices]

    # ---- Assemble feature matrix ----
    feat_list  = []
    feat_names = []
    if dist_arr is not None:
        feat_list.append(dist_arr)
        feat_names += [f"dist_{i}_{j}" for (i, j) in pair_indices]
    feat_list.append(rg_arr.reshape(-1,1)); feat_names.append("Rg")
    feat_list.append(rmsd_arr.reshape(-1,1)); feat_names.append("RMSD")
    feat_list.append(b_factors.reshape(-1,1)); feat_names.append("Bfactor")
    feat_list.append(avg_phi.reshape(-1,1)); feat_names.append("Avg_Phi")
    feat_list.append(avg_psi.reshape(-1,1)); feat_names.append("Avg_Psi")
    feat_list.append(frac_helix.reshape(-1,1)); feat_names.append("Frac_Helix")
    feat_list.append(frac_sheet.reshape(-1,1)); feat_names.append("Frac_Sheet")
    feat_list.append(frac_coil.reshape(-1,1)); feat_names.append("Frac_Coil")
    feat_list.append(ptm_contact_matrix); feat_names += [f"PTM_{lbl}" for lbl in ptm_site_labels]

    features = np.hstack(feat_list)
    print(f"Final feature shape: {features.shape}")
    return features, feat_names, ptm_contact_matrix, ptm_site_labels

def extract_and_save_features(
    system_name,
    info,
    step,
    feat_dir,
    radius_ptm,
    radius_nad,
    radius_interfaces,
    dssp_window,
    contact_cutoff
):
    u, frames = load_and_preprocess_subsystem_mda(system_name, info, step)
    sel = select_atoms_for_analysis(
        u,
        ptm_resnames      = info.get("ptm_resnames","AUTO"),
        radius_ptm        = radius_ptm,
        radius_nad        = radius_nad,
        radius_interfaces = radius_interfaces
    )
    if sel["final"].n_atoms == 0:
        feats      = np.zeros((len(frames), 1))
        names      = ["ZeroFeature"]
        ptm_mat    = np.zeros((len(frames), 1))
        ptm_labels = ["ZeroFeature"]
    else:
        feats, names, ptm_mat, ptm_labels = extract_features_mda(
            u,
            sel["final"],
            frames,
            compute_pairwise_dist=False,
            native_contact_cutoff=contact_cutoff,
            dssp_window=dssp_window,
            ptm_resnames_for_contacts=info.get("ptm_resnames",[])
        )
    os.makedirs(feat_dir, exist_ok=True)
    np.save(os.path.join(feat_dir, f"{system_name}_feats.npy"), feats)
    with open(os.path.join(feat_dir, f"{system_name}_featnames.json"), "w") as f:
        json.dump(names, f)
    # save PTM contacts & labels
    np.save(os.path.join(feat_dir, f"{system_name}_ptm_matrix.npy"), ptm_mat)
    with open(os.path.join(feat_dir, f"{system_name}_ptm_labels.json"), "w") as f:
        json.dump(ptm_labels, f)
    return system_name
