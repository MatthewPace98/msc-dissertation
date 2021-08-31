h=/home/mpace21/data

multiqc \
-o $h/multiqc_out \
$h/trim_galore_output/*.txt \
$h/fastqc/FastQC_out/* \
$h/fastqc/FastQC_filt_out/* \
$h/FastQ_Screen_out/filtered/* \
$h/FastQ_Screen_out/raw/*

