# Sample information
- Single-end reads
- Sanger/Illumina 1.9 encoding (ASCII+33 quality scores as Phred scores)
- 51bp reads
- 30x depth
- HiSeq 2000 Sequencing System (Illumina) used by EMBL in Heidelberg, Germany
- Total RNA was extracted from the treated cells using the RNeasy Mini kit (Qiagen) procedure

5 = Phenolic compound treated - 1 hr time point, 10 = Phenolic compound treated - 6 hr time point, 15 = Phenolic compound treated - 12 hr time point, 16 = Negative control untreated cells.

| File name                                                             | Treatment         | Size (GB) |
|-----------------------------------------------------------------------|-------------------|-----------|
| C9LMRACXX_5_16s004919-1-1_Zammit-Mangion_lane216s004919_sequence.txt  | 1 hour treatment  | 1.93      |
| C9LMRACXX_10_16s004924-1-1_Zammit-Mangion_lane316s004924_sequence.txt | 6 hour treatment  | 1.59      |
| C9LMRACXX_15_16s004929-1-1_Zammit-Mangion_lane416s004929_sequence.txt | 12 hour treatment | 1.30      |
| C9LMRACXX_16_16s004930-1-1_Zammit-Mangion_lane416s004930_sequence.txt | Negative control  | 1.48      |

# Preprocessing

## FastQC
http://manpages.ubuntu.com/manpages/impish/man1/fastqc.1.html

```fastqc -o [PATH] *.txt```

FastQC can accept multiple file names as input, so we can use the ```*.txt``` wildcard

Results were of typical RNA-seq data

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
If BWA not installed: ```sudo apt install bwa```                   

(in directory)  
```
perl fastq_screen \
--aligner bowtie2 \
--outdir /home/bioinf/Desktop/RNAseq/FASTQC/fastq_screen_output \
/home/bioinf/Desktop/RNAseq/FASTQ_files/C9LMRACXX_15_16s004929-1-1_Zammit-Mangion_lane416s004929_sequence.txt.gz
```

## PrinSeq++
                          

# Alignment - STAR

Genome index generation
```
STAR
```

Alignment
```
STAR
```
