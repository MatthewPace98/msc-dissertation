

## My DEGseq 

library("DEGseq")

home <- "/home/bioinf/Desktop/RNAseq"

one_bed <- file.path(home, "/STAR_out/1hourAligned.sortedByCoord.out.bed")
control_bed <- file.path(home, "/STAR_out/controlAligned.sortedByCoord.out.bed")
six_bed <- file.path(home, "/STAR_out/6hourAligned.sortedByCoord.out.bed")
twelve_bed <- file.path(home, "/STAR_out/12hourAligned.sortedByCoord.out.bed")

one_geneExp <- file.path(home, "RSEM_out/genes/1hour.genes.results")
control_geneExp <- file.path(home, "RSEM_out/genes/control.genes.results")
six_geneExp <- file.path(home, "RSEM_out/genes/6hour.genes.results")
twelve_geneExp <- file.path(home, "RSEM_out/genes/12hour.genes.results")

geneExpMatrix1 <- readGeneExp(file=one_geneExp, geneCol=1, valCol=c(6))
geneExpMatrix0 <- readGeneExp(file=control_geneExp, geneCol=1, valCol=c(6))
geneExpMatrix6 <- readGeneExp(file=six_geneExp, geneCol=1, valCol=c(6))
geneExpMatrix12 <- readGeneExp(file=twelve_geneExp, geneCol=1, valCol=c(6))

write.table(geneExpMatrix1[30:31,],row.names=FALSE)  #?
layout(matrix(c(1,2), 3, 2, byrow=TRUE))  #?
par(mar=c(2, 2, 2, 2))  #?

# Comparing control and 1hr
outputDir <- file.path(home, "DGEseq/prefilter", "1hr")
DEGexp(geneExpMatrix1=geneExpMatrix1, geneCol1=1, expCol1=c(2), groupLabel1="1hr", outputDir=outputDir, geneExpMatrix2=geneExpMatrix0, geneCol2=1, expCol2=c(2), groupLabel2="0hr", method="MARS")
# Comparing control and 6hr
outputDir <- file.path(home, "DGEseq/prefilter", "6hr")
DEGexp(geneExpMatrix1=geneExpMatrix6, geneCol1=1, expCol1=c(2), groupLabel1="6hr", outputDir=outputDir, geneExpMatrix2=geneExpMatrix0, geneCol2=1, expCol2=c(2), groupLabel2="0hr", method="MARS")
# Comparing control and 12hr
outputDir <- file.path(home, "DGEseq/prefilter", "12hr")
DEGexp(geneExpMatrix1=geneExpMatrix12, geneCol1=1, expCol1=c(2), groupLabel1="12hr", outputDir=outputDir, geneExpMatrix2=geneExpMatrix0, geneCol2=1, expCol2=c(2), groupLabel2="0hr", method="MARS")


# Comparing control and 1hr
outputDir <- file.path(home, "DGEseq/postfilter", "1hr")
DEGexp(geneExpMatrix1=geneExpMatrix1, geneCol1=1, expCol1=c(2), groupLabel1="1hr", outputDir=outputDir, geneExpMatrix2=geneExpMatrix0, geneCol2=1, expCol2=c(2), groupLabel2="0hr", method="MARS")
# Comparing control and 6hr
outputDir <- file.path(home, "DGEseq/postfilter", "6hr")
DEGexp(geneExpMatrix1=geneExpMatrix6, geneCol1=1, expCol1=c(2), groupLabel1="6hr", outputDir=outputDir, geneExpMatrix2=geneExpMatrix0, geneCol2=1, expCol2=c(2), groupLabel2="0hr", method="MARS")
# Comparing control and 12hr
outputDir <- file.path(home, "DGEseq/postfilter", "12hr")
DEGexp(geneExpMatrix1=geneExpMatrix12, geneCol1=1, expCol1=c(2), groupLabel1="12hr", outputDir=outputDir, geneExpMatrix2=geneExpMatrix0, geneCol2=1, expCol2=c(2), groupLabel2="0hr", method="MARS")

