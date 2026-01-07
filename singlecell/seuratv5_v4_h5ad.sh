> reticulate::py_config()
python:         /home/liuxiaodongLab/jiangjing/.cache/R/reticulate/uv/cache/archive-v0/_lFjfVE7RHx4dX9k2chwc/bin/python3
libpython:      /home/liuxiaodongLab/jiangjing/.cache/R/reticulate/uv/python/cpython-3.11.12-linux-x86_64-gnu/lib/libpython3.11.so
pythonhome:     /home/liuxiaodongLab/jiangjing/.cache/R/reticulate/uv/cache/archive-v0/_lFjfVE7RHx4dX9k2chwc:/home/liuxiaodongLab/jiangjing/.cache/R/reticulate/uv/cache/archive-v0/_lFjfVE7RHx4dX9k2chwc
virtualenv:     /home/liuxiaodongLab/jiangjing/.cache/R/reticulate/uv/cache/archive-v0/_lFjfVE7RHx4dX9k2chwc/bin/activate_this.py
version:        3.11.12 (main, Apr  9 2025, 04:04:00) [Clang 20.1.0 ]
numpy:          /home/liuxiaodongLab/jiangjing/.cache/R/reticulate/uv/cache/archive-v0/_lFjfVE7RHx4dX9k2chwc/lib/python3.11/site-packages/numpy
numpy_version:  2.2.5

NOTE: Python version was forced by py_require()
> reticulate::py_install("anndata", pip = TRUE)
> reticulate::py_module_available("anndata")
> library(sceasy)
> library(Seurat)
> data <- readRDS("/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/diffusionmap/new_umap_dim18_spread2.5_mindist0.7_celltypes_250416.rds")
data[["RNA"]] <- as(data[["RNA"]],"Assay")> data[["RNA"]] <- as(data[["RNA"]],"Assay")
> setwd("/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/diffusionmap")
> sceasy::convertFormat(data,from="seurat",to="anndata",outFile="new_umap_dim18_spread2.5_mindist0.7_celltypes_250416.h5ad")

