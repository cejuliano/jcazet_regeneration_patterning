#! /bin/bash -l
#SBATCH --job-name=consensus
#SBATCH -c 8
#SBATCH -t 7-0
#SBATCH --mem=50G
#SBATCH --error=consensus.err
#SBATCH --output=consensus.out

source resources/venv/bin/activate
module load samtools

resources/generate_consensus.sh
