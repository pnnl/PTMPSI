# plotting.py

import os
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib as mpl
import itertools
import matplotlib.gridspec as gridspec

from matplotlib.patches import FancyArrowPatch, Patch
from matplotlib.lines import Line2D
from adjustText import adjust_text

# import config to apply your global styling & access constants if needed
import pca_config


def plot_cluster_feature_ranges(cluster_feature_ranges, feat_names, system_name, outdir="figures_pca"):
    """
    Plots errorbar plots for the range (min and max) of each feature per cluster; Å->nm, Å²->nm²
    """
    units_dict = {
        "Rg": " Rg (nm)",
        "RMSD": " RMSD (nm)",
        "Bfactor": " B-factor (nm²)",
        "Avg_Phi": " Phi (rad)",
        "Avg_Psi": " Psi (rad)",
        "Frac_Helix": " Helix Count (%) ",
        "Frac_Sheet": " β-sheet Count (%) ",
        "Frac_Coil": " Coil Count (%) ",
    }
    
    n_feats = len(feat_names)
    n_cols = 3
    n_rows = int(np.ceil(n_feats / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    axes = axes.flatten()
    
    cluster_ids = sorted(cluster_feature_ranges.keys())
    for i, feat in enumerate(feat_names):
        ax = axes[i]
        means = []
        errors = []
        for cid in cluster_ids:
            if feat in cluster_feature_ranges[cid]:
                minv, maxv = cluster_feature_ranges[cid][feat]
                meanv = 0.5 * (minv + maxv)
                errv = 0.5 * (maxv - minv)
            else:
                meanv, errv = 0, 0
            means.append(meanv)
            errors.append(errv)
        xvals = np.arange(len(cluster_ids))
        # choose same colormap as plot_clusters_pca
        n_clusters = len(cluster_ids)
        cmap = plt.cm.get_cmap("tab10" if n_clusters <= 10 else "tab20", n_clusters)
        for j, cid in enumerate(cluster_ids):
            color = cmap(j)
            m, e = means[j], errors[j]
            ax.errorbar(xvals[j], m, yerr=e, fmt='o', capsize=3, color=color)
        ax.tick_params(labelsize=18)
        ax.set_xticks(xvals)
        ax.set_xticklabels([str(c) for c in cluster_ids])
        ax.set_title(feat, fontsize=18)
        unit_str = units_dict.get(feat, "")
        ax.set_ylabel(unit_str, fontsize=18)
        # ─── PTM‑contact features get their own xlabel ───
        if feat.startswith("PTM_"):
            label = feat[len("PTM_"):]
            resname = label.split("-")[0]
            ax.set_ylabel(f"Mean {resname} Residue Contacts")
            
    for j in range(i+1, len(axes)):
        axes[j].set_visible(False)
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.12)
    # add a single centered xlabel
    fig.text(
        0.5, 0.04, 
        "Cluster Index",
        ha="center", va="center",
        fontsize=mpl.rcParams["axes.labelsize"]
    )
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"{system_name}_cluster_feature_ranges.pdf")
    plt.savefig(outfile, dpi=300)
    print(f"Saved cluster feature ranges plot to {outfile}")
    plt.close('all')


def plot_top_pcs(pca_data, explained_variance, system_name, outdir="figures_pca"):

    top5 = min(pca_data.shape[1], 5)
    # use qualitative colormap 
    cmap = plt.cm.get_cmap("tab10" if top5 <= 10 else "tab20", top5)
    colors = [cmap(i) for i in range(top5)] # C0, C1, …

    fig = plt.figure(figsize=(20,8))
    gs = gridspec.GridSpec(2, 3, height_ratios=[2,1], width_ratios=[1,1,1])

    # Big PC1 panel
    ax_big = fig.add_subplot(gs[:,0])
    ax_big.hist(pca_data[:, 0], bins=30, alpha=0.7, color=colors[0])
    ax_big.set_title(f"PC1 ({explained_variance[0]*100:.1f}% var)", fontsize=18)
    ax_big.set_xlabel("PC1 Score", fontsize=18)
    ax_big.set_ylabel("Frequency", fontsize=18)

    # The smaller four PCs
    subaxes = [fig.add_subplot(gs[0,1]), fig.add_subplot(gs[0,2]),
               fig.add_subplot(gs[1,1]), fig.add_subplot(gs[1,2])]
    for i, ax in enumerate(subaxes, start=1):
        if i < top5:
            ax.hist(pca_data[:, i], bins=30, alpha=0.7, color=colors[i])
            ax.set_title(f"PC{i+1} ({explained_variance[i]*100:.1f}% var)", fontsize=18)
            ax.set_xlabel(f"PC{i+1} Score", fontsize=18)
            ax.set_ylabel("Frequency", fontsize=18)
        else:
            ax.set_visible(False)

    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"{system_name}_top5_PCs.pdf")
    
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    print(f"Saved top 5 PC components plot to {outfile}")
    plt.close('all')

def plot_clusters_pca(pca_data, cluster_labels, centroids, system_name, outdir="figures_pca", title="PCA Clusters"):
    n_clusters = len(np.unique(cluster_labels))
    # choose a qualitative colormap with exactly n_clusters colors
    cmap = plt.cm.get_cmap("tab10" if n_clusters <= 10 else "tab20", n_clusters)

    plt.figure(figsize=(8,6))
    # plot each cluster
    for c in range(n_clusters):
        idxs = np.where(cluster_labels == c)[0]
        plt.scatter(
            pca_data[idxs, 0], pca_data[idxs, 1],
            color=cmap(c),
            alpha=0.7,
            label=f"Cluster {c}"
        )

    # plot centroids
    for c, centroid in enumerate(centroids):
        plt.scatter(
            centroid[0], centroid[1],
            marker="X", s=200,
            color=cmap(c),
            edgecolor='black',
            linewidth=1.5
        )

    plt.xlabel("PC1", fontsize=18)
    plt.ylabel("PC2", fontsize=18)
    # legend
    plt.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, 1.10),
        ncol=min(n_clusters, 4),
        fontsize=18,
        frameon=False
    )

    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"{system_name}_PCA_clusters.pdf")
    plt.savefig(outfile, dpi=300)
    print(f"Saved PCA cluster plot to {outfile}")
    plt.close('all')

def visualize_ktn(G, T, mfpt, tau, system_name, dt_ns, cluster_labels=None, outdir="figures_pca", title="Kinetic Transition Network"):
    """
    Visualizes the KTN with node labels showing τ (in ns) and edge labels showing MFPT (in ns).
    Filters reciprocal edges so that each unordered pair is shown only once.
    """
    mfpt_ns = mfpt * dt_ns
    tau_ns  = tau * dt_ns
    
    # Filter edges: if both (i,j) and (j,i) exist, keep only the one for i<j
    filtered_edges = []
    for (i, j, d) in G.edges(data=True):
        if G.has_edge(j, i):
            if i < j:
                filtered_edges.append((i, j, d))
        else:
            filtered_edges.append((i, j, d))
    
    if cluster_labels is not None:
        uc, cts = np.unique(cluster_labels, return_counts=True)
        freq_dict = dict(zip(uc, cts))
        n_clusters = len(uc)
    else:
        freq_dict = {}
        n_clusters = G.number_of_nodes()

    # compute cluster percentages
    if cluster_labels is not None:
        total = len(cluster_labels)
        percentages = {node: freq_dict.get(node, 0) / total * 100 for node in G.nodes()}
    else:
        percentages = {node: 0.0 for node in G.nodes()}

    cmap = plt.cm.get_cmap("tab10" if n_clusters <= 10 else "tab20", n_clusters)
    node_colors = [cmap(node) for node in G.nodes()]
    
    pos = nx.spring_layout(G, seed=42)
    plt.figure(figsize=(14,14))
    
    # Draw edges with arrowheads exactly at node boundaries
    ax = plt.gca()
    # recompute node sizes as in your node-drawing call:
    sizes = [300 + 30*freq_dict.get(node,0) for node in G.nodes()]
    edge_arrows = []
    for (i, j, d) in filtered_edges:
        w = d['weight']
        r_src = np.sqrt(sizes[i])
        r_tgt = np.sqrt(sizes[j])
        arrow = FancyArrowPatch(
            pos[i], pos[j],
            arrowstyle='->',
            mutation_scale=30,
            shrinkA=0.5*r_src,
            shrinkB=0.5*r_tgt,
            linewidth=4.0*w,
            color='black',
            alpha=1.0,
            zorder=1
        )
        ax.add_patch(arrow)
        edge_arrows.append(arrow) 

    nodes_artist = nx.draw_networkx_nodes(
        G, pos,
        node_size=[300 + 30 * freq_dict.get(node, 0) for node in G.nodes()],
        node_color=node_colors,
        alpha=0.7,
        ax=ax
    )
    nodes_artist.set_zorder(2)
    
    node_text_objs = []
    edge_text_objs = []
    
    for node, (xx, yy) in pos.items():
        lbl_str = f"ID {node} ({percentages[node]:.1f}%)\nτ={tau_ns[node]:.2f} ns"
        t = plt.text(xx, yy, lbl_str, ha='center', va='top', color='black', zorder=3, fontsize=14)
        node_text_objs.append(t)
    
    for (i, j, d) in filtered_edges:
        prob = d['weight']
        xA, yA = pos[i]
        xB, yB = pos[j]
        xm = 0.2 * xB + 0.8 * xA
        ym = 0.2 * yB + 0.8 * yA
        lbl_str = f"p={prob:.2f}\nmfpt={mfpt_ns[i,j]:.2f} ns"
        t_edge = plt.text(xm, ym, lbl_str, color='darkblue', ha='center', va='top', zorder=3, fontsize=14)
        edge_text_objs.append(t_edge)
    
    adjust_text(
        node_text_objs, 
        arrowprops=dict(arrowstyle='-', color='gray', lw=0.5),
        only_move={'points':'xy','text':'xy','objects':'xy'}
    )
    adjust_text(
        node_text_objs + edge_text_objs,
        objects=edge_arrows,
        arrowprops=dict(arrowstyle='-', color='gray', lw=0.5),
        only_move={'points':'xy','text':'xy','objects':'xy'}
    )

    plt.axis('off')
    plt.margins(0.2)
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"{system_name}_KTN.pdf")
    plt.savefig(outfile, dpi=300, bbox_inches='tight', pad_inches=0.4)
    print(f"Saved KTN plot to {outfile}") 
    plt.close('all') 

def plot_cluster_frequency_3d(cluster_labels, overall_scores, system_name, outdir="figures_pca"):
    uc, counts = np.unique(cluster_labels, return_counts=True)
    avg_scores = [np.mean(overall_scores[cluster_labels == cval]) for cval in uc]
    n_clusters = len(uc)
    cmap = plt.cm.get_cmap("tab10" if n_clusters <= 10 else "tab20", n_clusters)

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')
    ax.tick_params(labelsize=18)
    for i, cval in enumerate(uc):
        ax.bar([cval], [counts[i]], zs=[avg_scores[i]], zdir='y', alpha=0.8, color=cmap(i))
    ax.set_xlabel("Cluster ID", fontsize=18)
    ax.set_ylabel("Avg Score", fontsize=18)
    ax.set_zlabel("Frequency", fontsize=18)
    outfile = os.path.join(outdir, f"{system_name}_cluster_frequency_3d.pdf")
    os.makedirs(outdir, exist_ok=True)
    plt.savefig(outfile, dpi=300)
    print(f"Saved 3D cluster freq plot to {outfile}")
    plt.close('all')

def plot_ptm_importance(ptm_site_labels, ptm_matrix, system_name, outdir="figures_pca"):
    """
    Plots the mean contacts for each PTM site.
    Also appends which domain the chainID belongs to (PRK, GAP, CP12) in the x-axis label.
    """
    # Helper: get domain from chain ID
    gap_segids  = ["A","D","E","F","I","K","N","O"]
    prk_segids  = ["B","G","J","L"]
    cp12_segids = ["C","H","M","P"]
    
    def domain_label(segid):
        if segid in gap_segids:
            return "GAP"
        elif segid in prk_segids:
            return "PRK"
        elif segid in cp12_segids:
            return "CP12"
        else:
            return "?"

    # —––––– compute means –––––—
    site_means = np.mean(ptm_matrix, axis=0)

    # —––––– assign tab10 colors in same order as wheel nodes, use tab20 if more than 10 sites, otherwise tab10 –––––—

    n_sites = len(ptm_site_labels)
    cmap_name = "tab10" if n_sites <= 10 else "tab20"
    cmap = plt.get_cmap(cmap_name, n_sites)
    bar_colors = [cmap(i) for i in range(n_sites)]

    # Build the x-tick labels with domain info
    labeled_sites = []
    for site_label in ptm_site_labels:
        parts = site_label.split("-")
        if len(parts) == 3:
            resname, segid, resid = parts
            labeled_sites.append(f"{resname}-{segid}-{resid} [{domain_label(segid)}]")
        else:
            labeled_sites.append(site_label)
    
    plt.figure(figsize=(8,5))
    plt.bar(labeled_sites, site_means, color=bar_colors, alpha=0.8)
    plt.xlabel("PTM Residue (resname-chainID-resid) [domain]", fontsize=18)
    plt.ylabel("Mean Residue Contacts", fontsize=18)
    plt.xticks(rotation=70, fontsize=18)
    outfile = os.path.join(outdir, f"{system_name}_ptm_importance.pdf")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    print(f"Saved PTM importance bar plot to {outfile}")
    plt.close('all')

def plot_top5_pcs_all_systems(pca_data_dict, outdir="figures_pca", filename="all_systems_top5_PCs.pdf"):
    systems = list(pca_data_dict.keys())
    # determine the maximum number of PCs across all systems
    max_pcs = max(pca_data.shape[1] for pca_data, _ in pca_data_dict.values())
    # choose a qualitative colormap with at least max_pcs colors
    cmap = plt.cm.get_cmap("tab10" if max_pcs <= 10 else "tab20", max_pcs)

    fig, axes = plt.subplots(len(systems), max_pcs, figsize=(4 * max_pcs, 3 * len(systems)), squeeze=False)
    for i, sys_name in enumerate(systems):
        pca_data, explained_variance = pca_data_dict[sys_name]
        n_pcs = pca_data.shape[1]

        # plot each available PC
        for j in range(n_pcs):
            ax = axes[i, j]
            ax.hist(pca_data[:, j], bins=30, alpha=0.7, color=cmap(j))
            ax.set_title(f"{sys_name} PC{j+1}\n({explained_variance[j]*100:.1f}% var)", fontsize=18)
            ax.set_xlabel(f"PC{j+1} Score", fontsize=18)
            ax.set_ylabel("Frequency", fontsize=18)

        # hide any unused subplots
        for j in range(n_pcs, max_pcs):
            axes[i, j].set_visible(False)

    plt.tight_layout()
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, filename)
    plt.savefig(outfile, dpi=300)
    print(f"Saved combined top 5 PCs plot for all systems to {outfile}")
    plt.close('all')

def plot_pairwise_top5_pcs(pca_data_dict, cluster_labels_dict=None, outdir="figures_pca", filename="pairwise_top5_PCs.pdf"):
    """
    For each system, create pairwise scatter plots of the top 5 PCs, 
    color-coded by cluster ID. Only add a legend to the first subplot,
    since the color code is the same for all subplots.
    """
    # Domain color not relevant here; we color by cluster ID
    for sys_name, (pca_data, explained_variance) in pca_data_dict.items():
        # If no cluster labels given, default to 0
        if cluster_labels_dict is not None and sys_name in cluster_labels_dict:
            cl_labels = cluster_labels_dict[sys_name]
        else:
            cl_labels = np.zeros(pca_data.shape[0], dtype=int)
        
        unique_clusters = np.unique(cl_labels)
        # pick tab10 for up to 10 clusters, otherwise tab20
        ncl = len(unique_clusters)
        cmap_name = 'tab10' if ncl <= 10 else 'tab20'
        color_map = plt.cm.get_cmap(cmap_name, ncl)
        
        pairs = list(itertools.combinations(range(5), 2))
        n_pairs = len(pairs)
        nrows, ncols = 2, 5
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
        axes = axes.flatten()
        
        # store scatter artists from the first subplot to build a legend
        scatter_handles = []
        scatter_labels  = []
        
        for subplot_index, (i, j) in enumerate(pairs):
            ax = axes[subplot_index]
            for cluster_id in unique_clusters:
                # points in that cluster
                idxs = np.where(cl_labels == cluster_id)[0]
                sc = ax.scatter(
                    pca_data[idxs, i], 
                    pca_data[idxs, j], 
                    alpha=0.7, 
                    color=color_map(cluster_id),
                    label=f"Cluster {cluster_id}" if subplot_index == 0 else None
                )
                # Only gather handles/labels in the first subplot
                if subplot_index == 0:
                    scatter_handles.append(sc)
                    scatter_labels.append(f"Cluster {cluster_id}")
            ax.set_xlabel(f"PC{i+1} ({explained_variance[i]*100:.1f}% var)", fontsize=18)
            ax.set_ylabel(f"PC{j+1} ({explained_variance[j]*100:.1f}% var)", fontsize=18)
            ax.set_title(f"{sys_name}: PC{i+1} vs PC{j+1}", fontsize=18)
        
        # Hide any extra axes (if pairs < nrows*ncols)
        for ax in axes[n_pairs:]:
            ax.set_visible(False)
        
        # Create a legend on the first subplot only
        # (or as a single legend for the figure)
        handles = scatter_handles  
        labels  = scatter_labels  
        ncol = len(handles)  
        fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.10), ncol=ncol,frameon=False, fontsize=18)   
        os.makedirs(outdir, exist_ok=True)
        outfile = os.path.join(outdir, filename)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(outfile, dpi=300)
        print(f"Saved pairwise top 5 PCs scatter plots for {sys_name} to {outfile}")
        plt.close('all')

def build_and_plot_residue_network(
    u,
    sel_dict,
    frames,
    system_name,
    group="ptm",
    cutoff=3.0,
    outdir="figures_pca"
):
    """
    PTM‐shell wheel network on a circle:
      - all PTM residues and their shell neighbors placed uniformly on a circle
      - edges connect each PTM residue to its neighbors
      - edge width ∝ average atom–atom contacts < cutoff Å over all frames
      - edges color‐coded by PTM origin
    """
    # pick the AtomGroup for this group
    if group == "ptm":
        total_ag = sel_dict.get("ptm_total", sel_dict["final"])
    elif group == "interface":
        total_ag = sel_dict["interfaces"]["prk_cp12"] + sel_dict["interfaces"]["gap_cp12"]
    elif group == "nad":
        total_ag = sel_dict.get("nad_total", sel_dict["final"])
    elif group == "nad_ptm":
        nad_sel = sel_dict.get("nad_total", u.atoms[[]])
        ptm_sel = sel_dict.get("ptm_total", u.atoms[[]])
        total_ag = nad_sel + ptm_sel
        if total_ag.n_atoms == 0:
            total_ag = sel_dict["final"]
    else:  # "all"
        total_ag = sel_dict["final"]
        for grp in sel_dict["interfaces"].values():
            total_ag += grp

    if total_ag.n_atoms == 0:
        print(f"No atoms in selection '{group}' for {system_name}; skipping network.")
        return

    # identify PTM residues (non‐standard)
    ptm_residues = [
        r for r in sel_dict["ptm_total"].residues
        if r.resname.upper() not in STANDARD_NAMES
    ]
    if not ptm_residues:
        print(f"No PTM residues for {system_name}; skipping network.")
        return

    # neighbors = residues in total_ag minus PTMs themselves
    ptm_keys = {(r.segid, r.resid) for r in ptm_residues}
    neighbors = [
        r for r in total_ag.residues
        if (r.segid, r.resid) not in ptm_keys
    ]
    if not neighbors:
        print(f"No PTM-shell neighbors for {system_name}; skipping network.")
        return

    # compute average contact counts using a single-array contact_matrix
    # flatten atom‐groups one time
    ptm_atoms = u.atoms[[atom.index for res in ptm_residues for atom in res.atoms]]
    nbr_atoms = u.atoms[[atom.index for res in neighbors     for atom in res.atoms]]

    # map each atom back to its PTM‐residue index
    ptm_res_idx_arr = np.repeat(np.arange(len(ptm_residues)),
                                [len(res.atoms) for res in ptm_residues])
    # map each atom back to its neighbor‐residue index
    nbr_res_idx_arr = np.repeat(np.arange(len(neighbors)),
                                [len(res.atoms) for res in neighbors])

    # initialize per‐residue contact counts
    counts = np.zeros((len(ptm_residues), len(neighbors)), dtype=float)

    for fi in frames:
        u.trajectory[fi]

        coords_ptm = ptm_atoms.positions
        coords_nbr = nbr_atoms.positions

        # stack them into one (N,3) array
        coords = np.vstack((coords_ptm, coords_nbr))
        n_ptm  = coords_ptm.shape[0]

        # calculate full NxN contact matrix
        # contact_matrix(coord, cutoff, returntype='numpy', box=...)
        cmat = contact_matrix(
            coords,
            cutoff,                    # here cutoff is the 2nd positional arg
            returntype="numpy",
            box=u.trajectory.ts.dimensions
        )

        # slice out only the PTM‐vs‐neighbor block
        mask = cmat[:n_ptm, n_ptm:]

        # if no contacts this frame, skip
        i_ptm, i_nbr = np.where(mask)
        if i_ptm.size == 0:
            continue

        # map atom‐indices back to residue‐indices
        res_i = ptm_res_idx_arr[i_ptm]
        res_j = nbr_res_idx_arr[i_nbr]

        # flatten residue‐pairs into a single bincount index
        flat_idx = res_i * counts.shape[1] + res_j
        additions = np.bincount(flat_idx, minlength=counts.size).reshape(counts.shape)
        counts += additions

    # normalize to average per frame
    counts /= len(frames)

    # rebuild edge_weights
    edge_weights = {}
    for i, ptm_res in enumerate(ptm_residues):
        clab = f"{ptm_res.resname}-{ptm_res.segid}-{ptm_res.resid}"
        for j, nbr_res in enumerate(neighbors):
            w = counts[i, j]
            if w > 0:
                nlab = f"{nbr_res.resname}-{nbr_res.segid}-{nbr_res.resid}"
                edge_weights[(clab, nlab)] = w

    if not edge_weights:
        print(f"No contacts found for {system_name}; skipping network.")
        return

    # build the graph
    G = nx.Graph()
    for center in ptm_residues:
        clab = f"{center.resname}-{center.segid}-{center.resid}"
        G.add_node(clab, type="ptm")
    for nbr in neighbors:
        nlab = f"{nbr.resname}-{nbr.segid}-{nbr.resid}"
        G.add_node(nlab, type="neighbor")
    for (clab, nlab), w in edge_weights.items():
        G.add_edge(clab, nlab, weight=w)

    # layout: circle
    all_nodes = [f"{r.resname}-{r.segid}-{r.resid}" for r in ptm_residues + neighbors]
    angles    = np.linspace(0, 2*np.pi, len(all_nodes), endpoint=False)
    pos       = {n:(np.cos(a),np.sin(a)) for n,a in zip(all_nodes, angles)}

    # draw
    plt.figure(figsize=(14,14))
    cmap = plt.get_cmap('tab10')
    ptm_nodes = [n for n,d in G.nodes(data=True) if d["type"]=="ptm"]
    nbr_nodes = [n for n,d in G.nodes(data=True) if d["type"]=="neighbor"]
    nx.draw_networkx_nodes(G,pos,nodelist=ptm_nodes,node_size=400,
                           node_color=[cmap(i) for i in range(len(ptm_nodes))],
                           label="PTM")
    nx.draw_networkx_nodes(G,pos,nodelist=nbr_nodes,node_size=200,
                           node_color="lightgray",label="Neighbor")

    max_w = max(edge_weights.values())
    ptm_color_map = {clab:cmap(i) for i,clab in enumerate(ptm_nodes)}
    for (clab,nlab),w in edge_weights.items():
        color = ptm_color_map.get(clab,'gray')
        width = 1 + 4*(w/max_w)
        nx.draw_networkx_edges(G,pos,edgelist=[(clab,nlab)],
                               width=width,edge_color=[color],alpha=0.7)

    # radial labels
    for node,(x,y) in pos.items():
        ang = np.degrees(np.arctan2(y,x))
        ha  = "left"
        if ang < -90 or ang > 90:
            ha = "right"; ang += 180
        plt.text(x,y,node,rotation=ang,rotation_mode="anchor",
                 horizontalalignment=ha,verticalalignment="center",
                 fontsize=18)

               
    plt.axis("off")

    # grab the existing PTM / neighbor handles
    node_handles, node_labels = plt.gca().get_legend_handles_labels()

    # create one colored‐line proxy per PTM residue
    edge_handles = [
        Line2D([0], [0], color=ptm_color_map[clab], lw=3)
        for clab in ptm_nodes
    ]
    edge_labels = ptm_nodes

    # combine node and edge entries into one legend
    all_handles = node_handles + edge_handles
    all_labels  = node_labels  + edge_labels
    ncol = math.ceil(len(all_handles)/2)

    plt.legend(
        all_handles,
        all_labels,
        loc='upper center',
        bbox_to_anchor=(0.5, 1.10),
        fontsize=18,
        ncol=ncol,
        frameon=False
    )

    # 8) save
    os.makedirs(outdir, exist_ok=True)
    outfn = os.path.join(outdir, f"{system_name}_ptm_wheel_network.pdf")
    plt.tight_layout()
    plt.savefig(outfn, dpi=300)
    plt.close('all')
    print("Saved PTM wheel network →", outfn)

