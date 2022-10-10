a=/home/mpace21/data
f=6hour
mkdir $a/prinseq_out

prinseq++ \
-fastq $a/trim_galore_output/${f}_trimmed.fq \
-out_name $a/prinseq_out/${f}_filtered.fq \
-out_bad $a/prinseq_out/${f}_bad.fq \
-ns_max_n 7 \
-lc_dust=0.1 > $a/prinseq_out/control_prinseq_report.txt

f=12hour

prinseq++ \
-fastq $a/trim_galore_output/${f}_trimmed.fq \
-out_name $a/prinseq_out/${f}_filtered.fq \
-out_bad $a/prinseq_out/${f}_bad.fq \
-ns_max_n 7 \
-lc_dust=0.1 > $a/prinseq_out/control_prinseq_report.txt

