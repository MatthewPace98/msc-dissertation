#!/usr/bin/env bash

# Generates a combined quality control report from Trim Galore!, FastQC, FastQScreen and RSEM output files.

h=/home/mpace21/data

multiqc \
-o $h/multiqc_out \
$h/trim_galore_output/*.txt \
$h/fastqc/FastQC_out/. \
$h/FastQ_Screen_out/raw/. \
$h/RSEM/out/.
