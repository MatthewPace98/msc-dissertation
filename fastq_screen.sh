#!/usr/bin/env bash

# fastq_screen aligns reads against a set of reference genomes to check data quality. 
# Input files are filtered reads from 'prinseq++' output.

a=/home/mpace21/data/FastQ_Screen_out

fastq_screen \
--aligner bowtie2 \
--outdir $a/filtered \
--conf $a/FastQ_Screen_Genomes/fastq_screen.conf \
/home/mpace21/data/prinseq_out/*good*

