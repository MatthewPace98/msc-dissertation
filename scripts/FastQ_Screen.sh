a=/home/mpace21/data/FastQ_Screen_out

fastq_screen \
--aligner bowtie2 \
--outdir $a/filtered \
--conf $a/FastQ_Screen_Genomes/fastq_screen.conf \
/home/mpace21/data/prinseq_out/*good*

