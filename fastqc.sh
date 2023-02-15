#!/usr/bin/env/ bash

# Analyses RNA sequence quality

h=/home/mpace21/data

fastqc \
-o $h/fastqc/FastQC_filt_out \
-t 8 \
/home/mpace21/data/prinseq_out/*good*
