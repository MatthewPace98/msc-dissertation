#!/bin/bash

# Performs quantification using RSEM on four fastq files.
# The tool is run with the STAR aligner and reference genome GRCh38, with the reads sorted by coordinate.
# Outputs gene expression estimates.

samples=(control 1hour 6hour 12hour)
input=/home/mpace21/data/prinseq_out/
h=/home/mpace21/data/RSEM

for sample in ${samples[@]};
do
	rsem-calculate-expression \
	--num-threads 16 \
	--star \
	--strandedness reverse \
	--sort-bam-by-coordinate \
	$input/${sample}*good_out.fastq \
	$h/reference/GRCh38 \
	$h/expression/${sample}

done
