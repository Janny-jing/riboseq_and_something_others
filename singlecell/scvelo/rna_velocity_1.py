#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import warnings
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend
import matplotlib.pyplot as plt
# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')
# Define JOB_ID
JOB_ID = "01"
# Cell 3: Define Paths
print("Defining paths...")

# Base paths
data_base_path = "/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/scvelo/data"
plots_base_path = "/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/scvelo/plot"

# Ensure the plots directory exists
os.makedirs(plots_base_path, exist_ok=True)

# File paths
velocyto_file = os.path.join(data_base_path, 'merged.loom')
clusters_file = os.path.join(data_base_path, "cell_clusters.csv")
umap_file = os.path.join(data_base_path, "cell_embeddings.csv")
path_10x = os.path.join(data_base_path, "filtered_feature_bc_matrix")

# Check if files exist
required_files = [clusters_file, umap_file, path_10x, velocyto_file]
for f in required_files:
    if not os.path.exists(f):
        print(f"Missing required file: {f}")
        raise FileNotFoundError(f"Missing data file: {f}")

print("All required files found.")
# Cell 4: Read Clusters
print("Reading clusters file...")
Clusters_Loupe = pd.read_csv(clusters_file, delimiter=',',index_col=0)
Barcodes = Clusters_Loupe.index
# Read UMAP exported from Loupe Browser 
UMAP_Loupe = pd.read_csv(umap_file, delimiter=',',index_col=0)
# Tansform to Numpy (for formatting)
UMAP_Loupe = UMAP_Loupe.to_numpy()
Sample3p = sc.read_10x_mtx(path_10x, var_names='gene_symbols')
Sample3p_df = Sample3p.to_df()
Sample3p = Sample3p[Barcodes]
# Add Clusters from Loupe to object
Sample3p.obs['Loupe'] = Clusters_Loupe

# Add UMAP from Loupe to object
Sample3p.obsm["X_umap"] = UMAP_Loupe
# Read velocyto output
VelNeutro3p = sc.read(velocyto_file)
VelNeutro3p.var = VelNeutro3p.var.set_index("var_names")
VelNeutro3p.obs = VelNeutro3p.obs.set_index("obs_names")

# Step 2: 清理索引名称（可选）
VelNeutro3p.var_names.name = None
VelNeutro3p.obs_names.name = None

# Step 3: 处理重复项（如有需要）
VelNeutro3p.var_names_make_unique()
VelNeutro3p.obs_names_make_unique()
Sample3p_union = scv.utils.merge(Sample3p, VelNeutro3p)
print("Computing velocities...")
scv.pp.filter_and_normalize(Sample3p_union, min_shared_counts=10, n_top_genes=500)
scv.pp.moments(Sample3p_union, n_pcs=30, n_neighbors=30)
#scv.tl.recover_dynamics(Sample3p_union, n_jobs=-1)
scv.tl.velocity(Sample3p_union, mode='steady_state')
scv.tl.velocity_graph(Sample3p_union, n_jobs=-1)
#scv.tl.recover_latent_time(Sample3p_union)
print("Velocities computed.")
Sample3p_union.write("velocity_results_1.h5ad")

