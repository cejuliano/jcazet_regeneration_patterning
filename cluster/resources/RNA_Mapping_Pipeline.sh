#! /bin/bash

# This code takes raw RNA-seq reads in a fastq.gz format and generates
# an RSEM results file, containing aligned read counts per gene. This
# output is then used 

# This code requires:

# Trimmomatic (V. 0.36)
# rsem (V. 1.2.31)
# Bowtie2 (V. 2.2.6)


# Read files should be formated as such:
# prefix_RNA.fastq.gz

# Prefix is the unique sample identifier that has the following format:
# Time[Treatment]StructureRep

# As an example, the first biological replicate for head regenerates collected 
# at 12 hours post amputation treated with iCRT14 would have the prefix:
# 12iH1

# The untreated equivalent would have the prefix:
# 12H1

########################

# The only argument needed for this script is the unique sample ID
prefix="$1"

echo "$prefix"

# create an initial fastqc report of both forward and reverse reads
# we'll compare this to the filtered reads to ensure that trimming 
# was effective 
echo "initial fastqc report"
fastqc -o . resources/"$prefix"_RNA*.fastq.gz

# Some of the samples have paired end data. In these cases we only 
# want to use read 1 for downstream analyses

# use trimmomatic to remove adapter sequences and
# stretches of low quality base calls.
# It also sets the minimum read size to 32
echo "filtering reads"

# Some of the samples have paired end reads. In these cases we will 
# only use read 1 for downstream analyses

# Some other samples were sequenced over multiple lanes
# For these we'll do the filtering and post filtering QC
# separately, then pool them for mapping

if [ -f resources/"$prefix"_RNA_R1.fastq.gz ]; then

    java -jar resources/trimmomatic-0.36.jar SE -threads 16  -phred33 \
        resources/"$prefix"_RNA_R1.fastq.gz "$prefix"_RNA_trim.fastq.gz \
        ILLUMINACLIP:resources/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 MINLEN:36
        
elif [ -f resources/"$prefix"_RNA_L3* ]; then

	for arg in resources/"$prefix"_RNA_L*
	do
	
		echo "$arg"
		outputName="${arg/.fastq.gz/_trim.fastq.gz}"
		outputName="${outputName/resources\//}"
		
		java -jar resources/trimmomatic-0.36.jar SE -threads 16  -phred33 \
       		"$arg" "$outputName" \
        	ILLUMINACLIP:resources/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 \
        	SLIDINGWINDOW:4:15 MINLEN:36
        	
    done
    
else

		java -jar resources/trimmomatic-0.36.jar SE -threads 16  -phred33 \
        resources/"$prefix"_RNA.fastq.gz "$prefix"_RNA_trim.fastq.gz \
        ILLUMINACLIP:resources/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 MINLEN:36
        
fi
     
# create a fastqc report of filtered reads
# so we can make sure filtering worked well

# you should see: high phred scores across the entire
# read length and no overrepresented sequences
echo "post-filtering fastqc report" 
fastqc -o . "$prefix"*trim.fastq.gz

# If reads are split over multiple lanes, pool them here

if [ -f "$prefix"_RNA_L3*trim.fastq.gz ]; then

	echo "pooling lanes"
	zcat "$prefix"_RNA_L*trim.fastq.gz | gzip > "$prefix"_RNA_trim.fastq.gz
	rm "$prefix"_RNA_L*trim.fastq.gz
	
fi

# We then map reads to the genome and count the number of
# aligned reads per transcript

echo "mapping reads"

# Rsem will output a genes.results file that will be used downstream to
# generate the read count matrix
rsem-calculate-expression --num-threads 16 --no-bam-output --bowtie2 \
	--temporary-folder "$prefix"_tmp "$prefix"_RNA_trim.fastq.gz \
        resources/Dovetail_mRNAs "$prefix"_RNA
        
rm "$prefix"*isoforms.results
rm -r "$prefix"*.stat

echo "done"
