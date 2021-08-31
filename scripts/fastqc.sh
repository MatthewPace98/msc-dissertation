#!/usr/bin/env/ bash

h=/home/mpace21/data

fastqc \
-o $h/fastqc/FastQC_filt_out \
-t 8 \
/home/mpace21/data/prinseq_out/*good*

# -a $h/FastQC/truseq_adapters.txt \
