# Preprocessing

## FastQC
http://manpages.ubuntu.com/manpages/impish/man1/fastqc.1.html

```
fastqc \
-o ./FastQC_out \
-t 6 \
*.fq.gz
```

FastQC can accept multiple file names as input, so we can use the ```.fq.gz``` wildcard

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

If reference genomes absent: ```fastq_screen --get_genomes```
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
https://github.com/Adrian-Cantu/PRINSEQ-plus-plus

```
prinseq++ \
-fastq $a/trim_galore_output/${f}_trimmed.fq \
-out_name $a/prinseq_out/${f}_filtered.fq \
-out_bad $a/prinseq_out/${f}_bad.fq \
-ns_max_n 7 \
-lc_dust=0.1 > $a/prinseq_out/control_prinseq_report.txt
```


## MultiQC
Dependencies:
- Python version 2.7+, 3.4+ or 3.5+

https://multiqc.info/docs/
In directory:
```multiqc .```


# Alignment and Quantification - STAR and RSEM
The trimmed FASTQ files were aligned to the GRCh38.p13 reference genome.

##Reference generation
```
rsem-prepare-reference --gtf mm9.gtf \
                       --star \
                       --star-sjdboverhang 50 \
                       --num-threads 8 \
                       --gtf $h2/Homo_sapiens.GRCh38.104.gtf \
                       ./Homo_sapiens.GRCh38.dna.primary_assembly.fa \
                       ./RSEM/reference/GRCh38
```

##Alignment
```
samples=(control, 1hour, 6hour, 12hour)
input=./prinseq_out/

for sample in ${samples[@]};
do
        rsem-calculate-expression \
        --num-threads 16 \
        --star \
        --sort-bam-by-coordinate \
        --fragment-length-mean 50 \
        --fragment-length-sd 1 \
        $input/${sample}*good_out.fq \
        ./reference/GRCh38 \
        ./expression/${sample}

done
```
