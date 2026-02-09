#!/bin/bash
#SBATCH --job-name bcdag_600
#SBATCH -q primary
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ht2699@wayne.edu
#SBATCH -o output_%j.out
#SBATCH -e errors_%j.err
#SBATCH -t 1-0:0:0

module load gnu7
module swap gnu7 gnu9
module load cmake/3.21.1
module load r/4.5.0

Rscript bcdag_ge_N_600_thinning.R


