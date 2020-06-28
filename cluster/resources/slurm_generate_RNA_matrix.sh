#!/bin/bash -l
#SBATCH --job-name=matrix
#SBATCH -c 1
#SBATCH -t 4-0
#SBATCH --mem=8G
#SBATCH --error=matrix.err
#SBATCH --output=matrix.out

module load rsem/1.2.31

resources/generate_RNA_Matrix.sh
