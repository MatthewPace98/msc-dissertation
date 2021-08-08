trim_galore \
--phred33 \
--fastqc \
-a "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" \
-o /home/bioinf/Desktop/RNAseq/trim_galore_output \
--fastqc_args "--outdir /home/bioinf/Desktop/RNAseq/FASTQC/trim_galore_output" \
/home/bioinf/Desktop/RNAseq/FASTQ_files/C9LMRACXX*

-q/--quality <INT>      Trim low-quality ends from reads in addition to adapter removal. For
                        RRBS samples, quality trimming will be performed first, and adapter
                        trimming is carried in a second round. Other files are quality and adapter
                        trimmed in a single pass. The algorithm is the same as the one used by BWA
                        (Subtract INT from all qualities; compute partial sums from all indices
                        to the end of the sequence; cut sequence at the index at which the sum is
                        minimal). Default Phred score: 20.

--phred33               Instructs Cutadapt to use ASCII+33 quality scores as Phred scores
                        (Sanger/Illumina 1.9+ encoding) for quality trimming. Default: ON.

--fastqc_args "<ARGS>"  Passes extra arguments to FastQC. If more than one argument is to be passed
                        to FastQC they must be in the form "arg1 arg2 etc.". An example would be:
                        --fastqc_args "--nogroup --outdir /home/". Passing extra arguments will
                        automatically invoke FastQC, so --fastqc does not have to be specified
                        separately.


-a/--adapter <STRING>   Adapter sequence to be trimmed. If not specified explicitly, Trim Galore will
                        try to auto-detect whether the Illumina universal, Nextera transposase or Illumina
                        small RNA adapter sequence was used. Also see '--illumina', '--nextera' and
                        '--small_rna'. If no adapter can be detected within the first 1 million sequences
                        of the first file specified or if there is a tie between several adapter sequences,
                        Trim Galore defaults to '--illumina' (as long as the Illumina adapter was one ofthe
                        options, else '--nextera' is the default). A single base
                        may also be given as e.g. -a A{10}, to be expanded to -a AAAAAAAAAA.

-o/--output_dir <DIR>   If specified all output will be written to this directory instead of the current
                        directory. If the directory doesn't exist it will be created for you.
                        

--length <INT>          Discard reads that became shorter than length INT because of either
                        quality or adapter trimming. A value of '0' effectively disables
                        this behaviour. Default: 20 bp.

                        For paired-end files, both reads of a read-pair need to be longer than
                        <INT> bp to be printed out to validated paired-end files (see option --paired).
                        If only one read became too short there is the possibility of keeping such                                                
                        unpaired single-end reads (see --retain_unpaired). Default pair-cutoff: 20 bp.


##########################################################
###################### FASTQ SCREEN ######################
##########################################################

# If BWA not installed:
sudo apt install bwa                        

(in directory)                        
perl fastq_screen \
--aligner bowtie2 \
--outdir /home/bioinf/Desktop/RNAseq/FASTQC/fastq_screen_output \
/home/bioinf/Desktop/RNAseq/FASTQ_files/C9LMRACXX_15_16s004929-1-1_Zammit-Mangion_lane416s004929_sequence.txt.gz

