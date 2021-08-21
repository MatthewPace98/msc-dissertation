# Preprocessing

## FastQC
http://manpages.ubuntu.com/manpages/impish/man1/fastqc.1.html

```fastqc -o ./FastQC_out *txt```

FastQC can accept multiple file names as input, so we can use the ```*.txt``` wildcard

Results were of typical RNA-seq data
FastQC was detecting TruSeq adapters 
Adapter info: 
http://tucf-genomics.tufts.edu/documents/protocols/TUCF_Understanding_Illumina_TruSeq_Adapters.pdf
https://dnatech.genomecenter.ucdavis.edu/wp-content/uploads/2013/06/illumina-adapter-sequences_1000000002694-00.pdf

## Trim Galore!
https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md

```
trim_galore \
--phred33 \
--fastqc \
-a "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" \
-o /home/bioinf/Desktop/RNAseq/trim_galore_output \
--fastqc_args "--outdir /home/bioinf/Desktop/RNAseq/FASTQC/trim_galore_output" \
/home/bioinf/Desktop/RNAseq/FASTQ_files/C9LMRACXX*
```

## FastQ Screen

https://github.com/StevenWingett/FastQ-Screen/blob/master/fastq_screen_documentation.md

Dependencies:
- perl
- BWA


(in directory)  
```
perl fastq_screen \
--aligner bowtie2 \
--outdir /home/bioinf/Desktop/RNAseq/FASTQC/fastq_screen_output \
/home/bioinf/Desktop/RNAseq/FASTQ_files/C9LMRACXX_15_16s004929-1-1_Zammit-Mangion_lane416s004929_sequence.txt.gz
```

## PrinSeq++
                          

## MultiQC
Dependencies:
- Python version 2.7+, 3.4+ or 3.5+

https://multiqc.info/docs/
In directory:
```multiqc .```


# Alignment - STAR

Genome index generation
```
STAR
```

Alignment
```
STAR
```
