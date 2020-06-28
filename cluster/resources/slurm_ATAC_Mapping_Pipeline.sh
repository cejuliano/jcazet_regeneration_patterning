#!/bin/bash -l
#SBATCH --job-name=ATAC_Pipeline
#SBATCH -c 16
#SBATCH -t 7-0
#SBATCH --array=0-70
#SBATCH --mem=50G
#SBATCH --error=ATAC_Pipeline_%a.err
#SBATCH --output=ATAC_Pipeline_%a.out

module load bowtie2
module load samtools
module load fastqc/0.11.4
module load bedtools
module load java/1.8

array=(0F1 \
0F2 \
0F3 \
0F4 \
0F5 \
0H1 \
0H2 \
0H3 \
0H4 \
0H5 \
0iF1 \
0iF2 \
0iF3 \
0iF4 \
0iH1 \
0iH2 \
0iH3 \
0iH4 \
12F1 \
12F2 \
12F3 \
12F4 \
12F5 \
12H1 \
12H2 \
12H3 \
12H4 \
12H5 \
12iF1 \
12iF2 \
12iF3 \
12iF4 \
12iF5 \
12iH1 \
12iH2 \
12iH3 \
12iH4 \
3F1 \
3F2 \
3F3 \
3H1 \
3H2 \
3H3 \
3iF1 \
3iF2 \
3iF3 \
3iF4 \
3iF5 \
3iH1 \
3iH2 \
3iH3 \
3iH4 \
8F1 \
8F2 \
8F3 \
8F4 \
8F5 \
8H1 \
8H2 \
8H3 \
8H4 \
8H5 \
8iF1 \
8iF2 \
8iF3 \
8iF4 \
8iH1 \
8iH2 \
8iH3 \
8iH4 \
8iH5 \
)

./resources/ATAC_Mapping_Pipeline.sh ${array[$SLURM_ARRAY_TASK_ID]}
