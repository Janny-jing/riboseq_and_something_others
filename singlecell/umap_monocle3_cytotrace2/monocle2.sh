#!/bin/bash
#SBATCH -c 2
#SBATCH -J cyto
#SBATCH --mem=64G
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal

module load R/4.4.3 gcc/11.2.0
Rscript difgene_monocle2_ncenter_dimension.R 
