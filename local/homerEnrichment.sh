#! /bin/bash

source ~/.bash_profile

bed2pos.pl "$1" > ePeak.peakfile.txt
bed2pos.pl "$2" > cPeak.peakfile.txt

findMotifsGenome.pl ePeak.peakfile.txt resources/Hm105_Dovetail_Assembly_1.0.fasta \
	Analysis_Output/ATAC/"$3" -size given -bg cPeak.peakfile.txt -nomotif -bits \
	-mknown resources/chromVar_HOMER.motifs -p 6
	
rm *peakfile.txt