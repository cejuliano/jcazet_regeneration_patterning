#! /bin/bash -l

shopt -s nullglob

# The only argument needed for this script is the unique sample ID
prefix="$1"

#index bam files
echo "indexing bam files"
for arg in "$prefix"[^_]*_final.bam;
do
	echo "$arg"
	samtools index "$arg"
done

echo "done"

#namesort bam files
echo "shifting reads"

for arg in "$prefix"[^_]*_final.bam;
do
	rep="${arg/_final.bam/}"
	echo "$rep"

	alignmentSieve -b "$arg" -o "$rep"_final_shift.bam -p 16 --ATACshift

	echo "resorting shifted reads"

	samtools sort -T "$rep".sort -o "$rep"_final_shift.sort.bam "$rep"_final_shift.bam

	mv "$rep"_final_shift.sort.bam "$rep"_final_shift.bam

	samtools index "$rep"_final_shift.bam
done

echo "done"

#make self pseudoreplicates
echo "making self-pseudoreplicates"

for arg in "$prefix"[^_]*_final_shift.bam;
do
	rep="${arg/_final_shift.bam/}"
	echo "$rep"

	echo "generating first psuedoreplicate"
	samtools view -s 1234.5 -b -@ 16 -o "$rep"_PR1_final_shift.bam "$arg"

	echo "generating second psuedoreplicate"
	samtools view "$rep"_PR1_final_shift.bam | cut -f 1 > "$rep"_PR1_qname.txt

	java -Xmx32g -jar resources/picard.jar FilterSamReads \
		I="$arg" \
		O="$rep"_PR2_final_shift.bam \
		READ_LIST_FILE="$rep"_PR1_qname.txt \
		VALIDATION_STRINGENCY=SILENT \
		FILTER=excludeReadList \
		QUIET=true

	rm "$rep"_PR1_qname.txt

done

#pool replicates

echo "pooling replicates"

samtools merge "$prefix"_MG_final_shift.bam "$prefix"[1-9]_final_shift.bam
samtools index "$prefix"_MG_final_shift.bam


#split pooled rep into psuedoreps (same number as total reps)

echo "generating pseudoreplicates from pooled counts"

numReps=("$prefix"[^_]*_final.bam)
numReps=${#numReps[@]}

echo "splitting into $numReps files"

count=1

cp "$prefix"_MG_final_shift.bam "$prefix"_MG_final_shift.sub.bam

while [ $numReps -gt 1 ]
do

	subSampleValue="$(Rscript resources/generateSubsampleValue.R $numReps | cut -f 2)"
	echo "$numReps"
	echo "$count"
	echo "$subSampleValue"
	samtools view -s "$subSampleValue" -b -@ 16 \
		-o "$prefix"_MG_PR"$count"_final_shift.bam "$prefix"_MG_final_shift.sub.bam

	samtools view "$prefix"_MG_PR"$count"_final_shift.bam | cut -f 1 > "$prefix"_PR_qname.txt


	java -Xmx32g -jar resources/picard.jar FilterSamReads \
		I= "$prefix"_MG_final_shift.sub.bam \
		O= "$prefix"_MG_final_shift.sub.tmp.bam \
		READ_LIST_FILE="$prefix"_PR_qname.txt \
		VALIDATION_STRINGENCY=SILENT \
		FILTER=excludeReadList \
		QUIET=true

	rm "$prefix"_MG_final_shift.sub.bam

	mv "$prefix"_MG_final_shift.sub.tmp.bam "$prefix"_MG_final_shift.sub.bam

	count=$(( $count + 1 ))
	numReps=$(( $numReps - 1 ))

done

mv "$prefix"_MG_final_shift.sub.bam "$prefix"_MG_PR"$count"_final_shift.bam

rm "$prefix"_PR_qname.txt

echo "Calling peaks"

for arg in "$prefix"*_final_shift.bam;
do
	rep="${arg/_final_shift.bam/}"
	echo "$rep"
	macs2 callpeak \
		-t "$arg" -f BAMPE -n "$rep" -g 9e8 -p 0.1 \
		--nomodel --keep-dup all

	sort -k 8gr,8gr "$rep"_peaks.narrowPeak \
		| awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' > "$rep".narrowPeak

	rm -f "$rep"_peaks.narrowPeak "$rep"_peaks.xls

done

rm *summits*

echo "done"

#perform idr for all reps
echo "Performing IDR on biological replicates"

for i in "$prefix"[1-9].narrowPeak
do
	for j in "$prefix"[1-9].narrowPeak
	do
		if [[ "$i" < "$j" ]]; then
			echo "$i"
			echo "$j"

			inputFile1="$i"
			prefix1="${inputFile1/.narrowPeak/}"

			inputFile2="$j"
			prefix2="${inputFile2/.narrowPeak/}"

			idr --samples "$i" "$j" \
				--peak-list "$prefix"_MG.narrowPeak --input-file-type narrowPeak \
				--output-file "$prefix1"_"$prefix2".idr --rank p.value --soft-idr-threshold 0.1

			IDR_THRESH_TRANSFORMED=$(awk -v p=0.1 'BEGIN{print -log(p)/log(10)}')

			awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' \
				{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' "$prefix1"_"$prefix2".idr | \
					sort | \
					uniq | \
					sort -k7n,7n > "$prefix1"_"$prefix2".IDR.narrowPeak
		fi
	done
done

echo "done"

#perform idr for all self-pseudoreps
echo "Perfoming IDR on self-pseudoreplicates"

for arg in "$prefix"[1-9]_final_shift.bam
do
	rep="${arg/_final_shift.bam/}"
	echo "$rep"

	inputFile1="$rep"_PR1.narrowPeak
	prefix1="${inputFile1/.narrowPeak/}"

	inputFile2="$rep"_PR2.narrowPeak
	prefix2="${inputFile2/.narrowPeak/}"

	idr --samples "$inputFile1" "$inputFile2" \
		--peak-list "$rep".narrowPeak --input-file-type narrowPeak \
		--output-file "$prefix1"_"$prefix2".idr --rank p.value --soft-idr-threshold 0.1

	IDR_THRESH_TRANSFORMED=$(awk -v p=0.1 'BEGIN{print -log(p)/log(10)}')

	awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' \
		{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' "$prefix1"_"$prefix2".idr | \
			sort | \
			uniq | \
			sort -k7n,7n > "$prefix1"_"$prefix2".IDR.narrowPeak
done

echo "done"

#perform idr for all pooled pseudoreps
echo "Perfoming IDR on pooled pseudoreplicates"

for i in "$prefix"_MG_PR*.narrowPeak
do
	for j in "$prefix"_MG_PR*.narrowPeak
	do
		if [[ "$i" < "$j" ]]; then
			echo "$i"
			echo "$j"

			inputFile1="$i"
			prefix1="${inputFile1/.narrowPeak/}"

			inputFile2="$j"
			prefix2="${inputFile2/.narrowPeak/}"

			idr --samples "$i" "$j" \
				--peak-list "$prefix"_MG.narrowPeak --input-file-type narrowPeak \
				--output-file "$prefix1"_"$prefix2".idr --rank p.value --soft-idr-threshold 0.1

			IDR_THRESH_TRANSFORMED=$(awk -v p=0.1 'BEGIN{print -log(p)/log(10)}')

			awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' \
				{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' "$prefix1"_"$prefix2".idr | \
					sort | \
					uniq | \
					sort -k7n,7n > "$prefix1"_"$prefix2".IDR.narrowPeak
		fi
	done
done

echo "done"

#use diffbind to generate consensus peakset
echo "Finding consensus peaksets"

echo "loading peaklists"
Rscript resources/Consensus_Peaks.R "$prefix"

echo "done"

#calculate reproducibility stats
echo "Calculating IDR rescue statistics"


echo "Rescue Ratio" > "$prefix"_IDR_Stats.txt

psLength=$(wc -l "$prefix"_psReps.consensus.bed | cut -f1 -d ' ')
bioRLength=$(wc -l "$prefix"_bioReps.consensus.bed | cut -f1 -d ' ')

if [[ $psLength > $bioRLength ]]
then

	echo "Ideal Replication Peakset Length/True Replication Peakset Length" \
        	>> "$prefix"_IDR_Stats.txt
	
	echo $psLength "/" $bioRLength >> "$prefix"_IDR_Stats.txt

	rescueRatio=$(echo "scale=3;$psLength/$bioRLength" | bc)

else

	echo "True Replication Peakset Length/Ideal Replication Peakset Length" \
		>> "$prefix"_IDR_Stats.txt

	echo $bioRLength "/" $psLength

	rescueRatio=$(echo "scale=3;$bioRLength/$psLength" | bc)
fi

echo "$rescueRatio" >> "$prefix"_IDR_Stats.txt

echo "" >> "$prefix"_IDR_Stats.txt

echo "Self-Consistency Ratio" >> "$prefix"_IDR_Stats.txt

echo "repA ideal peakset length" "/" "repB ideal peakset length" >> "$prefix"_IDR_Stats.txt

for i in "$prefix"[1-9]_PR*IDR.narrowPeak
do
	for j in "$prefix"[1-9]_PR*IDR.narrowPeak
	do
		iLength="$(wc -l $i | cut -f1 -d ' ')"
		jLength="$(wc -l $j | cut -f1 -d ' ')"
		if [[ "$i" < "$j" ]]; then
			if [[ "$iLength" < "$jLength" ]]; then

				echo "$j" "/" "$i" >> "$prefix"_IDR_Stats.txt
				echo "$jLength" "/" "$iLength" >> "$prefix"_IDR_Stats.txt
				selfRatio=$(echo "scale=3;$jLength/$iLength" | bc)
				echo "$selfRatio" >> "$prefix"_IDR_Stats.txt
				echo "" >> "$prefix"_IDR_Stats.txt

			else

				echo "$i" "/" "$j" >> "$prefix"_IDR_Stats.txt
				echo "$iLength" "/" "$jLength" >> "$prefix"_IDR_Stats.txt
				selfRatio=$(echo "scale=3;$iLength/$jLength" | bc)
				echo "$selfRatio" >> "$prefix"_IDR_Stats.txt
				echo "" >> "$prefix"_IDR_Stats.txt

			fi
		fi
	done
done

echo "done"

#generate bigwig files
echo "Generating bigwig tracks"

for arg in "$prefix"[1-9]_final_shift.bam
do
	echo "$arg"
	bamCoverage -b "$arg" -o "${arg/.bam/.bw}" \
		-of "bigwig" -bs 10 -p 16 --normalizeUsing "CPM"
done

echo "$prefix"_MG_final_shift.bam
bamCoverage -b "$prefix"_MG_final_shift.bam -o "$prefix"_MG_final_shift.bw \
	-of "bigwig" -bs 10 -p 16 --normalizeUsing "CPM"

echo "done"

#calculate TSS for all replicates

echo "Calculating TSS enrichment scores"

for arg in "$prefix"[1-9]_final_shift.bw
do
	rep="${arg/_final_shift.bw/}"

	computeMatrix reference-point --referencePoint TSS \
		-b 1000 -a 1000 \
		-R resources/dovetail_genes_High2000.bed \
		-S "$arg" -o "$rep"_ATAC_TSS_matrix.gz \
		--outFileSortedRegions regions_"$rep"_ATAC.bed \
		--outFileNameMatrix Values_"$rep"_ATAC.txt \
		--missingDataAsZero -p 16
done

Rscript resources/TSS_Calculation.R "$prefix"

echo "done"

echo "cleaning up"

rm regions_"$prefix"*ATAC.bed Values_"$prefix"*.txt "$prefix"*matrix.gz "$prefix"*.idr

echo "done"

