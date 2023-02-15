#!/bin/bash

# Prepares a reference genome index for use with rsem-calculate-expression 

h=/home/mpace21/data
h2=$h/miscfiles

rsem-prepare-reference --gtf mm9.gtf \
                       --star \
                       --star-sjdboverhang 50 \
                       --num-threads 8 \
                       --gtf $h2/Homo_sapiens.GRCh38.104.gtf \
                       $h2/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
                       $h/RSEM/reference/GRCh38
