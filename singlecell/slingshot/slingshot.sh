#!/bin/bash
#SBATCH -c 1
#SBATCH -J slingshot
#SBATCH --mem=64G
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal

module load R/4.4.3
Rscript slingshot.R 
