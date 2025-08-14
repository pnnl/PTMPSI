# pca_clustering.py

import numpy as np
import networkx as nx

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA, IncrementalPCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


def perform_pca(features, n_components):
    # ——— drop dummy all‐zero feature if present ———
    if features.shape[1] == 1 and np.all(features == 0):
        zeros = np.zeros((features.shape[0], 1))
        pca_model = PCA(n_components=1, random_state=42)
        pca_model.fit(zeros)
        return zeros, pca_model

    scaler = StandardScaler()
    scaled = scaler.fit_transform(features)
    pca_model = PCA(n_components=n_components, random_state=42)
    pca_data = pca_model.fit_transform(scaled)
    print("PCA done. shape:", pca_data.shape)
    print("Explained variance ratio:", pca_model.explained_variance_ratio_)
    return pca_data, pca_model

def cluster_kmeans(pca_data, n_clusters):
    km = KMeans(n_clusters=n_clusters, n_init="auto", random_state=42)
    labels = km.fit_predict(pca_data)
    print(f"Clustering K-Means into {n_clusters} clusters... Done.")
    return labels, km

def build_transition_matrix(cluster_labels, n_clusters):
    counts = np.zeros((n_clusters, n_clusters), dtype=int)
    for i in range(len(cluster_labels)-1):
        counts[cluster_labels[i], cluster_labels[i+1]] += 1
    T = counts.astype(float)
    for i in range(n_clusters):
        s = T[i].sum()
        if s > 0:
            T[i] /= s
    return T

def build_ktn(T):
    G = nx.DiGraph()
    n = T.shape[0]
    for i in range(n):
        G.add_node(i)
    for i in range(n):
        for j in range(n):
            if T[i, j] > 0:
                G.add_edge(i, j, weight=T[i, j])
    return G

def compute_stationary_distribution(T):
    eigvals, eigvecs = np.linalg.eig(T.T)
    idx = np.argmin(np.abs(eigvals - 1))
    pi = np.real(eigvecs[:, idx])
    pi /= np.sum(pi)
    return pi

def compute_cluster_score(T):
    return compute_stationary_distribution(T)

def compute_mfpt(T):
    """
    Mean first‐passage time matrix for transition matrix T.
    MFPT[i,j] = expected time to hit state j starting from i.
    If stationary probability pi[j] == 0, MFPT[i,j] is set to np.inf.
    """
    # 1) stationary distribution
    pi = compute_stationary_distribution(T)

    # 2) fundamental matrix Z = (I - T + E)^(-1), with E = 1·piᵀ
    n = T.shape[0]
    I = np.eye(n)
    E = np.outer(np.ones(n), pi)
    Z = np.linalg.inv(I - T + E)

    # 3) build MFPT matrix
    mfpt = np.zeros_like(T, dtype=float)
    for i in range(n):
        for j in range(n):
            if i != j:
                if pi[j] > 0:
                    mfpt[i, j] = (Z[j, j] - Z[i, j]) / pi[j]
                else:
                    mfpt[i, j] = np.inf
    return mfpt

def compute_residence_times(T):
    tau = np.zeros(T.shape[0])
    for i in range(T.shape[0]):
        if T[i,i] < 1:
            tau[i] = 1 / (1 - T[i,i])
        else:
            tau[i] = float('inf')
    return tau

def compute_local_score(pca_data, cluster_labels, kmeans_model, ptm_matrix=None):
    """
    Example local score: a combination of the sum of absolute PC coords + distance to cluster centroid,
    optionally minus some small penalty from PTM contacts (if provided).
    """
    sum_abs_pc = np.sum(np.abs(pca_data), axis=1)
    centroids = kmeans_model.cluster_centers_
    dists = np.array([np.linalg.norm(vec - centroids[label]) for label, vec in zip(cluster_labels, pca_data)])
    base_score = 0.5 * sum_abs_pc + 0.5 * dists
    if ptm_matrix is None:
        return base_score
    else:
        return base_score - 0.1 * np.sum(ptm_matrix, axis=1)

def overall_score(local_score, cluster_score, cluster_labels, alpha=0.6, beta=0.4):
    """
    Weighted combination of local_score and cluster_score (stationary distribution).
    """
    res = np.zeros_like(local_score)
    for i, ls in enumerate(local_score):
        res[i] = alpha * ls + beta * cluster_score[cluster_labels[i]]
    return res

def compute_combined_score(
    pc1_values: np.ndarray,
    cluster_labels: np.ndarray,
    mfpt: np.ndarray,
    taus: np.ndarray,
    global_scaler: StandardScaler,
    alpha: float = 1,
    beta:  float = 1,
    gamma: float = 1
) -> np.ndarray:
    """
    Combined per‐frame score = α·(z‐scored PC1) 
                            + β·(z‐scored mean‐outgoing‐MFPT) 
                            + γ·(z‐scored residence‐time).
    """
    # Build raw arrays
    mfpt_by_frame = np.array([np.nanmean(mfpt[c, :]) for c in cluster_labels])
    tau_by_frame  = np.array([taus[c]                for c in cluster_labels])

    # stack and clamp infinities/NaNs
    X = np.vstack((pc1_values, mfpt_by_frame, tau_by_frame)).T  # shape (n_frames,3)
    # replace inf/−inf with the largest/smallest finite value, and NaN with 0
    finite_mask = np.isfinite(X)
    if not finite_mask.all():
        max_finite = np.max(X[finite_mask])
        min_finite = np.min(X[finite_mask])
        X[~finite_mask & (X > 0)] = max_finite
        X[~finite_mask & (X < 0)] = min_finite
    X = np.nan_to_num(X, nan=0.0)

    # z‐score using the global scaler
    Xz = global_scaler.transform(X)
    pc1_z, mfpt_z, tau_z = Xz.T

    # Return weighted sum
    return alpha * pc1_z + beta * mfpt_z + gamma * tau_z
