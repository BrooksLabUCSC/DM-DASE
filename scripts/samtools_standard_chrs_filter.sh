#!/bin/bash

cd /PATH_TO_DATA/Aligned_BAM_Files

export PATH="PATH_TO_ENV/miniconda2/bin:$PATH"
export PATH="PATH_TO_ENV/bin/bedtools-2.25.0/bin:$PATH"
source PATH_TO_ENV/bin/activate

# Filter standard chrs only
for index in sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
	do
	echo "Processing "${index}".sortedByCoord.out.bam ..."
	samtools view -b "${index}".sortedByCoord.out.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > "${index}".sortedByCoord.trimmed.bam
	echo "Processing "${index}".sortedByCoord.out.bam completed..."
done

