#! /bin/bash

blastp -num_threads 6 -outfmt 6 -max_target_seqs 2 \
	-db nve -out hvToNve.txt -query hvCandidates.fa

blastp -num_threads 6 -outfmt 6 -max_target_seqs 2 \
	-db humanSP -out hvToHs.txt -query hvCandidates.fa
	
blastp -num_threads 6 -outfmt 6 -max_target_seqs 2 \
	-db epa -out hvToEpa.txt -query hvCandidates.fa