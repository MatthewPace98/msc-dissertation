#!/usr/bin/env/ bash

h=/home/mpace21/data/STAR

t=32 # Number of threads
h1=/home/mpace21/data/STAR
sample=1hour_filtered.fq

echo 'Aligning...'

STAR \
--runThreadN $t \
--genomeLoad NoSharedMemory \
--genomeDir $h1/index \
--readFilesIn $s1/home/mpace21/data/prinseq_out/${sample}_filtered_good_out.fastq \
--outFileNamePrefix $h1/out/${sample}_ \
--outSAMattributes NH HI AS nM NM MD RG \
--outFilterIntronMotifs RemoveNoncanonical \
--outSAMattrIHstart 0 \
--outSAMtype BAM Unsorted \
--outSAMattrRGline 'ID':${sample} \
--sjdbGTFfile $h1/Homo_sapiens.GRCh38.104.gtf \
--quantMode GeneCounts

echo 'Done'

# --genomeLoad >>> controls how the genome is loaded into memory. Options:
# - LoadAndKeep: load genome into shared and keep it in memory after run.
# - LoadAndRemove: load genome into shared but remove it after run.
# - LoadAndExit: load genome into shared memory and exit, keeping the genome in memory for future runs.
# - Remove: do not map anything, just remove loaded genome from memory.
# - NoSharedMemory: do not use shared memory, each job will have its own private copy of the genome.
#
# --readFilesIn [/path/to/read1] [/path/to/read2 ] >>> files containing the sequences to be mapped
#
# If the read files are compressed, use the --readFilesCommand UncompressionCommand option, where UncompressionCommand is the
# un-compression command that takes the file name as input parameter, and sends the uncompressed output to stdout. For example,
# for gzipped files (*.gz) use --readFilesCommand zcat OR --readFilesCommand gunzip -c. For bzip2- compressed files, use
# --readFilesCommand bunzip2 -c.
#
# For paired-end reads, use comma separated list for read1, followed by space, followed by comma separated list for read2, e.g.:
# --readFilesIn s1read1.fq,s2read1.fq,s3read1.fq s1read2.fq,s2read2.fq,s3read2.fq
# For multiple read files, the corresponding read groups can be supplied with space/comma/space- separated list in
# --outSAMattrRGline, e.g. --outSAMattrRGline ID:sample1 , ID:sample2 , ID:sample3
# Note that this list is separated by commas surrounded by spaces (unlike --readFilesIn list).
#
# --outFileNamePrefix >>> STAR produces multiple output files. All files have standard name, however, you can change the file
# prefixes using --outFileNamePrefix /path/to/output/dir/prefix.
#
# --outSAMattributes NH ...  >>> A string of desired SAM attributes, in the order desired for the output SAM.
# Tags can be listed in any combination/order.
# NH: number of loci the reads maps to: =1 for unique mappers, >1 for multimappers.
# HI: multiple alignment index, starts with –outSAMattrIHstart (=1 by default).
# AS: local alignment score, +1/-1 for matches/mismatches, score* penalties for indels and gaps. For PE reads, total score for two>
# nM: number of mismatches. For PE reads, sum over two mates.
# NM: edit distance to the reference (number of mismatched + inserted + deleted bases) for each mate.
# MD: string encoding mismatched and deleted reference bases (see standard SAM specifications).
# RG: read group. The read group to which the read belongs. If @RG headers are present, then readgroup must match the RG-ID
# field of one of the headers.
#
# --outFilterIntronMotifs RemoveNoncanonical >>> filter out alignments that contain non-canonical junctions.
# It is recommended to remove the non-canonical junctions for Cufflinks runs.
#
# --outSAMattrIHstart >>> The number of loci Nmap a read maps to is given by NH:i:Nmap field. Value of 1 corresponds to unique map>
#



