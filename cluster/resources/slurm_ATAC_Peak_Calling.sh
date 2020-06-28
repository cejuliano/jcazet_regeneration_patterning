#! /bin/bash -l
#SBATCH --job-name=ATAC_Peaks
#SBATCH -c 16
#SBATCH -t 7-0
#SBATCH --array=0-15
#SBATCH --mem=32G
#SBATCH --error=ATAC_Peaks_%a.err
#SBATCH --output=ATAC_Peaks_%a.out

source resources/venv/bin/activate
module load samtools

array=(0H 0F 3H 3F 8H 8F 12H 12F \
	0iH 0iF 3iH 3iF 8iH 8iF 12iH 12iF)


resources/ATAC_Peak_Pipeline.sh ${array[$SLURM_ARRAY_TASK_ID]}
