#!/usr/bin/env bash

h=/home/mpace21/data/STAR

mkdir $h/genomeDir

# Genome generation
STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir $h/index \
--genomeFastaFiles $h/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile $h/Homo_sapiens.GRCh38.104.gtf \
--sjdbOverhang 50

# runMode genomeGenerate >>> directs STAR to run genome indices generation job.
# --genomeFastaFiles >>>specifies genome reference sequence(s)
# --sjdbGTFfile >>>annotation file
# --sjdbOverhang  >>>read length minus 1
# --genomeDir  >>>where genome index will be generated

# --genomeDir >>> specifies path to the directory where the genome indices will be stored.
# This directory has to be created (with mkdir)before STAR runs and needs to have writing permissions.
# The file system needs to have at least 100GB of disk space available for a typical mammalian genome.
# It is recommended to remove all files from the genome directory before running the genome generation step.
# This directory path will have to be supplied at the mapping step to identify the reference genome.

#--genomeLoad >>> controls how the genome is loaded into memory.
# Options:
#- LoadAndKeep: load genome into shared and keep it in memory after run.
#- LoadAndRemove: load genome into shared but remove it after run.
#- LoadAndExit: load genome into shared memory and exit, keeping the genome in memory for future runs.
#- Remove: do not map anything, just remove loaded genome from memory.
#- NoSharedMemory: do not use shared memory, each job will have its own private copy of the genome.

#--genomeFastaFiles >>> specifies one or more FASTA files with the genome reference sequences.
#Multiple reference sequences (henceforth called “chromosomes”) are allowed for each fasta file.
#You can rename the chromosomes’ names in the chrName.txt keeping the order of the chromosomes in the file:
#the names from this file will be used in all output alignment files (such as.sam).
#The tabs are not allowed in chromosomes’ names, and spaces are not recommended.

#--sjdbGTFfile >>> specifies the path to the file with annotated transcripts in the standard GTF format.
# STAR will extract splice junctions from this file and use them to greatly improve accuracy of the mapping.
# While this is optional, and STAR can be run without annotations, using annotations is highly recommended whenever they are
# available.
#Starting from 2.4.1a, the annotations can also be included on the fly at the mapping step.

# --sjdbOverhang >>> specifies the length of the genomic sequence around the annotated junction to be used in constructing the
# splice junctions database. Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads.
# For instance, for Illumina 2x100b paired-end reads, the ideal value is 100-1=99.
# In case of reads of varying length, the ideal value is max(ReadLength)-1.
# In most cases, the default value of 100 will work as well as the ideal value.
