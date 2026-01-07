#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
from scipy.sparse import issparse
import os
import matplotlib.pyplot as plt
import numpy as np

def preprocess_data(file_path: str, resolution: float = 1.0, output_dir: str = None) -> dict:
    """
    Preprocesses single-cell data for downstream analysis.
    
    Parameters:
    -----------
    file_path : str
        Path to the input .h5ad file.
    resolution : float, optional
        Resolution parameter for clustering, controls the number of clusters.
        Default is 1.0.
    output_dir : str, optional
        Directory to save preprocessed outputs. If None, creates directory 
        based on input file location.
    
    Returns:
    --------
    dict
        Dictionary containing preprocessing results and summary information.
    """
    try:
        # Create output folder path
        file_name = os.path.basename(file_path)
        base_filename = os.path.splitext(file_name)[0]  # Remove .h5ad extension
        
        if output_dir is None or output_dir == "" or output_dir.strip() == "":
            # Use same logic as pages/preprocess_data.py: user's outputs folder
            try:
                # First try streamlit session state
                import streamlit as st
                if hasattr(st.session_state, 'user_id') and st.session_state.user_id:
                    user_id = st.session_state.user_id
                    output_dir = f"users/{user_id}/outputs"
                    print(f"Using user outputs directory from session: {output_dir}")
                else:
                    raise Exception("No user_id in session state")
            except:
                # Second try: look for any existing user directory (agent execution context)
                users_base = "users"
                if os.path.exists(users_base):
                    user_dirs = [d for d in os.listdir(users_base) if os.path.isdir(os.path.join(users_base, d))]
                    if user_dirs:
                        # Use the first user directory found
                        user_id = user_dirs[0]
                        output_dir = f"users/{user_id}/outputs"
                        print(f"Using detected user outputs directory: {output_dir}")
                    else:
                        # Final fallback - create a default user directory
                        output_dir = "users/guest/outputs"
                        print(f"Using default guest outputs directory: {output_dir}")
                else:
                    # Last resort fallback to _analysis (should rarely happen now)
                    parent_dir = os.path.dirname(file_path) 
                    output_dir = os.path.join(parent_dir, f"{base_filename}_analysis")
                    print(f"Final fallback to analysis directory: {output_dir}")
        else:
            print(f"Using provided output directory: {output_dir}")
        
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Create preprocessed output filename
        output_filename = f"{base_filename}_preprocessed.h5ad"
        
        # Load query data
        try:
            query_adata = sc.read_h5ad(file_path)
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {file_path}")
        
        print("Query dataset loaded successfully.")
        
        # Sanitize column names in .obs and .var to avoid reserved names
        if "_index" in query_adata.obs.columns:
            query_adata.obs.rename(columns={"_index": "cell_index"}, inplace=True)
        if "_index" in query_adata.var.columns:
            query_adata.var.rename(columns={"_index": "gene_index"}, inplace=True)
        

        # Remove "empty" genes
        sc.pp.filter_genes(query_adata, min_cells=1)

        # Convert adata.X to a numpy array if it's sparse or has 'todense'/'toarray' methods
        if hasattr(query_adata.X, 'todense'):
            X_data = np.asarray(query_adata.X.todense())  # Convert to dense matrix and then to NumPy array
        elif hasattr(query_adata.X, 'toarray'):
            X_data = np.asarray(query_adata.X.toarray())  # Convert to dense matrix and then to NumPy array
        else:
            X_data = np.asarray(query_adata.X)  # Already dense, convert to NumPy array

        # Check for NaN values
        if np.isnan(X_data).any():
            nan_count = np.isnan(X_data).sum()
            print(f"Found {nan_count} NaN values in data matrix")
            replace_nans = True  # Set to False if you want to skip instead

            if replace_nans:
                print(f"Replacing NaN values with zeros...")
                X_data = np.nan_to_num(X_data, nan=0.0)
                query_adata.X = X_data
            else:
                print(f"NaN values found in data, skipping...")
                return {
                    "status": "error",
                    "message": f"Found {nan_count} NaN values in data matrix",
                    "file_path": file_path
                }

        # Preprocess query dataset
        sc.settings.seed = 42  # Set random seed for reproducibility
        
        # Save raw counts
        query_adata.layers["counts"] = query_adata.X.copy()
        
        # Normalize total counts
        sc.pp.normalize_total(query_adata, target_sum=1e4)
        
        # Log-transform the data
        sc.pp.log1p(query_adata)
        query_adata.layers["logcounts"] = query_adata.X.copy()
        
        # Identify highly variable genes
        sc.pp.highly_variable_genes(query_adata, n_top_genes=2000, flavor="cell_ranger")
        
        # Perform PCA
        sc.tl.pca(query_adata, n_comps=30, use_highly_variable=True)
        
        # Compute neighborhood graph
        sc.pp.neighbors(query_adata)

        # Generate UMAP embeddings
        sc.tl.umap(query_adata)
        
        # Cluster cells using Leiden algorithm
        print(f"Clustering with resolution {resolution}...")
        sc.tl.leiden(query_adata, resolution=resolution, key_added=f"leiden_r{resolution}", random_state=42)
        
        
        # Export UMAP plot with dataset-specific name
        umap_plot_filename = f"{base_filename}_umap_plot_resolution_{resolution}.png"
        umap_plot_path = os.path.join(output_dir, umap_plot_filename)
        
        # Create UMAP plot
        plt.figure(figsize=(6, 4))
        sc.pl.umap(query_adata, color=[f"leiden_r{resolution}"], show=False, frameon=False)  
        plt.savefig(umap_plot_path, bbox_inches='tight')
        plt.close()
        print(f"UMAP plot saved at: {umap_plot_path}")
        
        # Fix any potential _index issues
        if query_adata.raw is not None:
            if '_index' in query_adata.raw.var.columns:
                query_adata.raw.var.rename(columns={'_index': 'index'}, inplace=True)
        else:
            if '_index' in query_adata.var.columns:
                query_adata.var.rename(columns={'_index': 'index'}, inplace=True)

        # Remove sparsity
        def remove_sparsity(adata):
            if issparse(adata.X):
                adata.X = adata.X.toarray()
            return adata

        query_adata = remove_sparsity(query_adata)
        
        # Save preprocessed data
        preprocessed_file_path = os.path.join(output_dir, output_filename)
        query_adata.write_h5ad(preprocessed_file_path)
        
        # Export metadata to CSV for download
        metadata_filename = f"{base_filename}_metadata.csv"
        metadata_path = os.path.join(output_dir, metadata_filename)
        query_adata.obs.to_csv(metadata_path)
        print(f"Metadata exported to: {metadata_path}")
        
        # Determine unique cluster count
        unique_clusters = len(query_adata.obs[f"leiden_r{resolution}"].unique())
        
        print(f"✅ Preprocessing completed successfully!")
        print(f"Data saved to: {preprocessed_file_path}")
        print(f"UMAP plot saved to: {umap_plot_path}")
        print(f"Dataset shape: {query_adata.shape}")
        print(f"Number of clusters (resolution {resolution}): {unique_clusters}")
        
        return {
            "status": "success",
            "message": f"Preprocessing completed successfully! Dataset shape: {query_adata.shape}, Clusters: {unique_clusters}",
            "file_path": file_path,
            "preprocessed_file": preprocessed_file_path,
            "umap_plot": umap_plot_path,
            "metadata_csv": metadata_path,
            "output_dir": output_dir,
            "dataset_shape": query_adata.shape,
            "cluster_count": unique_clusters,
            "resolution": resolution
        }
        
    except Exception as e:
        return {
            "status": "error",
            "message": str(e),
            "file_path": file_path
        }

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Preprocess single-cell data')
    parser.add_argument('file_path', help='Path to .h5ad file')
    parser.add_argument('--resolution', type=float, default=1.0, help='Clustering resolution')
    parser.add_argument('--output_dir', type=str, help='Output directory for results')
    
    args = parser.parse_args()
    
    result = preprocess_data(args.file_path, args.resolution, args.output_dir)
    print(f"Result: {result}")