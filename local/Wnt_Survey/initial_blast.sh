#! /bin/bash
blastp -num_threads 6 -outfmt 6 -max_target_seqs 2 \
	-db hvProt -out hsToHv.txt -query hsaWntGenes.fasta

blastp -num_threads 6 -outfmt 6 -max_target_seqs 2 \
	-db hvProt -out hmToHv.txt -query hmgWntGenes.fasta
	
blastp -num_threads 6 -outfmt 6 -max_target_seqs 2 \
	-db hvProt -out nvToHv.txt -query nveWntGenes.fasta
	
blastp -num_threads 6 -outfmt 6 -max_target_seqs 2 \
	-db hvProt -out epToHv.txt -query epaWntGenes.fasta