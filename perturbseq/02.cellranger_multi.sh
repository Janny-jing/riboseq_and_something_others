#! /bin/sh

/home/jj2024/01.soft/cellranger-8.0.1/cellranger multi \
	  --id=crispr_gex_analysis \
          --csv=multi_config.csv \
          --localcores=16 \
          --localmem=64
