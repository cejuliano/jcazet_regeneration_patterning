#! /bin/bash

blastn -num_threads 6 -outfmt 6 -max_target_seqs 5 \
	-db hvNucl -out lrToDv.txt -query aepLRv2.fasta
	
blastn -num_threads 6 -outfmt 6 -max_target_seqs 5 \
	-db lr -out dvToLr.txt -query Dovetail_mRNAs.fa