# Preprocessing

## FastQC
http://manpages.ubuntu.com/manpages/impish/man1/fastqc.1.html

```fastqc \
-o ./FastQC_out \
-t 6 \
*.gz
```

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
--length 45 \
-o $a/trim_galore_output \
--fastqc_args "--outdir /home/bioinf/Desktop/RNAseq/trim_galore_output" \
$a/Raw/*.fq.gz
```

## FastQ Screen

https://github.com/StevenWingett/FastQ-Screen/blob/master/fastq_screen_documentation.md

If referene genomes absent: ```fastq_screen --get_genomes```
Config file was adjusted
Tested on recommended test dataset


Dependencies:
- perl
- Bowtie2

```
fastq_screen \
--aligner bowtie2 \
--outdir /home/bioinf/Desktop/RNAseq/fastq_screen_output \
/home/bioinf/Desktop/RNAseq/Raw/*fq.gz
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
