h=/home/mpace21/data

mkdir $h/trim_galore_output

trim_galore \
--phred33 \
--fastqc \
-a "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" \
--length 45 \
-o $h/trim_galore_output \
--fastqc_args "-o /home/mpace21/data/trim_galore_output -t 8" \
--cores 4 \
$h/FastQ_files/*.fq
