#! /bin/bash

# make untreated consensus
for arg in *[0-9][FH]_bioReps.consensus.bed
do
	echo "$arg"
done

cat *[0-9][FH]_bioReps.consensus.bed > untreated_consensus.bed.tmp

./resources/bedClip -truncate untreated_consensus.bed.tmp resources/Dovetail.genome untreated_consensus.clip.bed.tmp

sort -k1,1 -k2,2n untreated_consensus.clip.bed.tmp > untreated_consensus.clip.sort.bed.tmp

bedtools merge -i untreated_consensus.clip.sort.bed.tmp > untreated_consensus.tmp.bed

samtools merge -@ 6 -f untreated_MG_final_shift.bam *[0-9][FH]_MG_final_shift.bam

samtools index untreated_MG_final_shift.bam

macs2 callpeak -t untreated_MG_final_shift.bam -f BAMPE -n untreated_MG \
	-g 9e8 -p 0.1 --nomodel --keep-dup all

sort -k 8gr,8gr untreated_MG_peaks.narrowPeak \
	| awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' > untreated_MG.narrowPeak

bedtools intersect -u -f 0.5 -a untreated_MG.narrowPeak -b untreated_consensus.tmp.bed \
	> untreated_consensus.bed

rm -f untreated_MG_peaks.narrowPeak \
	untreated_MG_peaks.xls \
	untreated_MG.narrowPeak \
	untreated_consensus.tmp.bed

# make full consensus
for arg in *_bioReps.consensus.bed
do
        echo "$arg"
done

cat *_bioReps.consensus.bed > full_consensus.bed.tmp

./resources/bedClip -truncate full_consensus.bed.tmp resources/Dovetail.genome full_consensus.clip.bed.tmp

sort -k1,1 -k2,2n full_consensus.clip.bed.tmp > full_consensus.clip.sort.bed.tmp

bedtools merge -i full_consensus.clip.sort.bed.tmp > full_consensus.tmp.bed

samtools merge -@ 6 -f full_MG_final_shift.bam *[!d]_MG_final_shift.bam

samtools index full_MG_final_shift.bam

macs2 callpeak -t full_MG_final_shift.bam -f BAMPE -n full_MG \
        -g 9e8 -p 0.1 --nomodel --keep-dup all

sort -k 8gr,8gr full_MG_peaks.narrowPeak \
        | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' > full_MG.narrowPeak

bedtools intersect -u -f 0.5 -a full_MG.narrowPeak -b full_consensus.tmp.bed \
        > full_consensus.bed

rm -f full_MG_peaks.narrowPeak \
        full_MG_peaks.xls \
        full_MG.narrowPeak \
        full_consensus.tmp.bed

