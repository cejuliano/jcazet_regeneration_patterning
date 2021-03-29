#! /bin/bash

source ~/.bash_profile

for f in injuryClusts/*.motif; do

	echo "$f"
	outputFile="${f/.motif/.hits.bed}"
	peakHits="${outputFile/.bed/.peaks.bed}"
	
	scanMotifGenomeWide.pl "$f" ../Hm105_Dovetail_Assembly_1.0.fasta \
		-bed -p 6 > "$outputFile" 2>/dev/null
		
	bedtools intersect -wo -a "$outputFile" \
		-b ../untreated_consensus_diffbind_labels.bed > "$peakHits".tmp
	
	bedtools merge -c 4,10 -o distinct -i "$peakHits".tmp > "$peakHits"
	
	rm injuryClusts/*tmp
		
done