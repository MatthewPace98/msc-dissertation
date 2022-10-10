samples=(control 1hour 6hour 12hour)
input=/home/mpace21/data/prinseq_out/
h=/home/mpace21/data/RSEM

for sample in ${samples[@]};
do
	rsem-calculate-expression \
	--num-threads 16 \
	--star \
	--strandedness reverse \
	--sort-bam-by-coordinate \
	$input/${sample}*good_out.fastq \
	$h/reference/GRCh38 \
	$h/expression/${sample}

done



# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323
# https://deweylab.github.io/RSEM/README.html
#
# rsem-prepare-reference and rsem-calculate-expression automatically use RSEM as the aligner, --star will allow RSEM to use the STAR aligner
# rsem-calculate-expression - expression calculation
# --fragment-length-mean
# --fragment-length-sd
#  --alignments
#        Input file contains alignments in SAM/BAM/CRAM format. The exact file
#        format will be determined automatically. (Default: off)
# When --no-bam-output is not specified and --sort-bam-by-coordinate is specified, RSEM will produce these three files:
# sample_name.transcript.bam, the unsorted BAM file, sample_name.transcript.sorted.bam and sample_name.transcript.sorted.bam.bai
#
# A wiggle plot representing the expected number of reads overlapping each position in the genome/transcript set can be
# generated from the sorted genome/transcript BAM file output.
# rsem-bam2wig sorted_bam_input wig_output wiggle_name [--no-fractional-weight]
#
# rsem-plot-model sample_name output_plot_file
# sample_name: the name of the sample analyzed
# output_plot_file: the file name for plots generated from the model. It is a pdf file
#
