# config.py

import os
import glob
import io
import json
import time
import tempfile
import argparse
import cProfile
import pstats

import numpy as np
import pandas as pd

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import networkx as nx

import MDAnalysis as mda
from MDAnalysis import Universe, Writer
from MDAnalysis.coordinates.PDB import PDBWriter
from MDAnalysis.analysis import distances, align
from MDAnalysis.analysis.rms import rmsd as mda_rmsd

import mdtraj as md

from joblib import Parallel, delayed

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA, IncrementalPCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from tabulate import tabulate
from matplotlib.patches import FancyArrowPatch, Patch
from matplotlib.lines import Line2D


# ─── Global matplotlib styling ───────────────────────────────────────────────
mpl.rcParams.update({
    "font.size":      14,
    "font.family":    "sans-serif",
    "font.weight":    "bold",
    "legend.fontsize":14,
    "axes.titlesize": 14,
    "axes.labelsize": 14,
    "xtick.labelsize":14,
    "ytick.labelsize":14,
    "axes.labelpad":  12,
    "xtick.major.pad":6,
    "ytick.major.pad":6
})

# ─── Standard residue names (20 AA + variants + ligands/ions) ────────────────
STANDARD_NAMES = {
    "ABU","ACE","AIB","ALA","ARG","ARGN","ASN","ASN1","ASP","ASP1","ASPH","ASPP","ASH",
    "CT3","CYS","CYS1","CYS2","CYSH","DALA","GLN","GLU","GLUH","GLUP","GLH","GLY",
    "HIS","HIS1","HISA","HISB","HISH","HISD","HISE","HISP","HSD","HSE","HSP","HYP",
    "ILE","LEU","LSN","LYS","LYSH","MELEU","MET","MEVAL","NAC","NME","NHE","NH2","PHE",
    "PHEH","PHEU","PHL","PRO","SER","THR","TRP","TRPH","TRPU","TYR","TYRH","TYRU","VAL",
    "PGLU","HID","HIE","HIP","LYP","LYN","CYN","CYM","CYX","DAB","ORN","NALA","NGLY",
    "NSER","NTHR","NLEU","NILE","NVAL","NASN","NGLN","NARG","NHID","NHIE","NHIP",
    "NHISD","NHISE","NHISH","NTRP","NPHE","NTYR","NGLU","NASP","NLYS","NORN","NDAB",
    "NLYSN","NPRO","NHYP","NCYS","NCYS2","NMET","NASPH","NGLUH","CALA","CGLY","CSER",
    "CTHR","CLEU","CILE","CVAL","CASN","CGLN","CARG","CHID","CHIE","CHIP","CHISD",
    "CHISE","CHISH","CTRP","CPHE","CTYR","CGLU","CASP","CLYS","CORN","CDAB","CLYSN",
    "CPRO","CHYP","CCYS","CCYS2","CMET","CASPH","CGLUH",
    "HOH","WAT","TIP3","UNK","SOL",
    "NAD","NADH","NADPH","NAI",
    "MG","ZN","NA","CL","K","CA"
}
