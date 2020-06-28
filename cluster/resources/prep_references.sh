#! /bin/bash

cd resources

bowtie2-build 105_mitochondrial_genome.fa hydra_mito

bowtie2-build Hm105_Dovetail_Assembly_1.0.fasta hydra_genome

rsem-prepare-reference --bowtie2 --transcript-to-gene-map Dovetail_mRNAs_Genemap.txt Dovetail_mRNAs.fa Dovetail_mRNAs
