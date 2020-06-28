#!/bin/bash -l
#SBATCH --job-name=RNA_Pipeline
#SBATCH -c 16
#SBATCH -t 7-0
#SBATCH --array=0-71
#SBATCH --mem=50G
#SBATCH --error=RNA_Pipeline_%a.err
#SBATCH --output=RNA_Pipeline_%a.out

module load bowtie2
module load fastqc/0.11.4
module load rsem

array=(0F1 0F2 0F3 0H1 0H2 0H3 \
        3F1 3F2 3F3 3H1 3H2 3H3 \
        8F1 8F2 8F3 8H1 8H2 8H3 \
        12F1 12F2 12F3 12H1 12H2 12H3 \
        0iF1 0iF2 0iF3 0iH1 0iH2 0iH3 \
        8iF1 8iF2 8iF3 8iH1 8iH2 8iH3 \
        12iF1 12iF2 12iF3 12iH1 12iH2 12iH3 \
        0wF1 0wF2 0wF3 0wH1 0wH2 0wH3 \
        2wF1 2wF2 2wF3 2wH1 2wH2 2wH3 \
        4wF1 4wF2 4wF3 4wH1 4wH2 4wH3 \
        8wF1 8wF2 8wF3 8wH1 8wH2 8wH3 \
        16wF1 16wF2 16wF3 16wH1 16wH2 16wH3)

resources/RNA_Mapping_Pipeline.sh ${array[$SLURM_ARRAY_TASK_ID]}
