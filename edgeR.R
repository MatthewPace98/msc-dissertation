library(edgeR)
library(org.Hs.eg.db)
library(AnnotationDbi) 

# https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# https://www.bioconductor.org/packages/devel/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html

# Files need to contain two columns, one for the counts and one for a gene identifier.

files <- list("1hour.genes.results", "12hour.genes.results", "6hour.genes.results", "control.genes.results")
path <- "/home/bioinf/Desktop/RNAseq/RSEM_out/genes/"

# Creating a DGEList which can be manipulated like any other list in R
# Takes normalised (TPM) counts as input
# The main components of an DGEList object are a matrix counts containing the integer counts,
# a data.frame samples containing information about the samples or libraries, and a optional
# data.frame genes containing annotation for the genes or genomic features

labels = c('1hr', '12hr', '6hr', '0hr')
group <- factor(labels)
y <- readDGE(files, path=path, columns=c(1,6), group=group, labels=labels)


# Filters genes with low counts 
# Threshold of 5-10 is normally taken
keep <- filterByExpr(y, group=group)
table(keep)
y <- y[keep,keep.lib.sizes=FALSE]

# https://www.youtube.com/watch?v=8qvGNTVz3Ik
# input 'keys' i.e. genes.
# Specify the 'keytype', so ifthey are Entrez IDS, Ensembl etc
# Specify the output, so what you would like to convert them to
Symbol <- mapIds(org.Hs.eg.db, keys=rownames(y), keytype="ENSEMBL",
                 column=c("ENTREZID"))
y$genes <- data.frame(Symbol=Symbol)
head(y$genes)



# Normalisation by TMM
y <- calcNormFactors(y)
design <- model.matrix(~group)

# Produces a plot in which distances between samples
# correspond to leading biological coefficient of variation (BCV) between those samples
#png(file="/home/bioinf/Desktop/RNAseq/edgeR/MDS.png")
#plotMDS(y)
#dev.off()

# Dispersion cannot be calculated so we take the BCV as 0.1 since we have no replicates
bcv <- 0.1
et_12 <- exactTest(y, pair=c(1, 2), dispersion=bcv^2)
et_1 <- exactTest(y, pair=c(1, 3), dispersion=bcv^2)
et_6 <- exactTest(y, pair=c(1, 4), dispersion=bcv^2)

#Display top differentially expressed tags
# logFC, the log-abundance ratio, i.e. fold change, for each tag in the two groups being compared
# logCPM, the log-average concentration/abundance for each tag in the two groups being compared
# PValue, exact p-value for differential expression using the NB model
# FDR, the p-value adjusted for multiple testing as found using p.adjust using the method specified
# FDR is calculated using the method by Benjamini & Hochberg (1995)
# Note: This may mess up annotations
# p.value is filtering by adjusted p, not raw p (i.e. FDR)
et_1_toptags <- topTags(et_1, n=nrow(et_1$table), adjust.method="BH", sort.by="PValue", p.value=0.05)$table
et_12_toptags <- topTags(et_12, n=nrow(et_12$table), adjust.method="BH", sort.by="PValue", p.value=0.05)$table
et_6_toptags <- topTags(et_6, n=nrow(et_6$table), adjust.method="BH", sort.by="PValue", p.value=0.05)$table

# Replaces et with et_toptags
et_1$table <- et_1_toptags[2:4]
et_1$genes <- et_1_toptags[1]

et_12$table <- et_12_toptags[2:4]
et_12$genes <- et_12_toptags[1]

et_6$table <- et_6_toptags[2:4]
et_6$genes <- et_6_toptags[1]


# To perform quasi-likelihood F-tests:
#QLfit <- glmQLFit(y,design, dispersion=bcv^2)
#qlf <- glmQLFTest(fit,coef=4,  dispersion=bcv^2)
#qlf_toptags <- topTags(qlf, n=20, adjust.method="BH", sort.by="logFC")$table

# To perform likelihood ratio tests:
#fit <- glmFit(y,design, dispersion=bcv^2)
#lrt <- glmLRT(fit,coef=4)
#lrt_toptags <- topTags(lrt, n=1000, adjust.method="BH", sort.by="PValue")$table

# Plot log-fold change against log-counts per million, with DE genes highlighted
#png(file="/home/bioinf/Desktop/RNAseq/edgeR/MD_4.png")
#plotMD(et)
#dev.off()


# A commonly used approach is to conduct DE tests, apply a fold-change cut-off 
# and then rank all the genes above that fold-change threshold by p-value
# The following is a rigorous statistical test for thresholded hypotheses.
# it tests whether the log2-fold-change is greater than lfc in absolute value.
#fit <- glmQLFit(y, design, dispersion=bcv)
#tr <- glmTreat(fit, coef=ncol(fit), lfc=2)
#tr_toptags <- topTags(tr)$table


# KEGG and GO


geneid_12 <-et_12_toptags$Symbol
geneid_1 <-et_1_toptags$Symbol
geneid_6 <-et_6_toptags$Symbol

head(geneid_12)
head(geneid_1)
head(geneid_6)

go_12 <- goana(et_12, species="Hs", geneid=geneid_12)
keg_12 <- kegga(et_12, species="Hs", geneid=geneid_12)

go_1 <- goana(et_1, species="Hs", geneid=geneid_1)
keg_1 <- kegga(et_1, species="Hs", geneid=geneid_1)

go_6 <- goana(et_6, species="Hs", geneid=geneid_6)
keg_6 <- kegga(et_6, species="Hs", geneid=geneid_6)

# write.csv(y$counts,"/home/bioinf/Desktop/RNAseq/Top DE genes/filtered_counts.csv", row.names = TRUE)
# write.csv(lrt_toptags,"/home/bioinf/Desktop/RNAseq/Top DE genes/lrttop.csv", row.names = TRUE)
write.csv(et_toptags,"/home/bioinf/Desktop/RNAseq/Top DE genes/ettop6hr.csv", row.names = TRUE)
write.csv(go,"/home/bioinf/Desktop/RNAseq/Top DE genes/go.csv", row.names = TRUE)
write.csv(keg,"/home/bioinf/Desktop/RNAseq/Top DE genes/kegg.csv", row.names = TRUE)


library(pathview)
# https://pathview.r-forge.r-project.org/pathview.pdf
# This package can be divided into four functional modules: the Downloader, 
# Parser, Mapper and Viewer
# gene.data can hold one sample or multiple sampleswith gene IDs as row names
# pathway.id the KEGG pathway ID usually 5 digits
# species is either the kegg code or name of target species. 'ko' is used for 
# KEGG ortholog pathway kegg.native whether to render as png (TRUE) or graphiz 
# layout (FASLE)
# https://www.genome.jp/kegg/pathway.html

# check code backup to see wtf is going on
?pathview

path_fun <- function(et, n, pathwayid){
  genedata <- as.matrix(et$table[, 1])
  rownames(genedata) <- et$genes[, 1]
  head(genedata)
  '
  https://www.genome.jp/kegg/pathway.html#metabolism
  05221 Acute myeloid leukemia
  05200 Pathways in cancer
  05202 Transcriptional misregulation in cancer
  04210 Apoptosis
  04110 Cell cycle
  04640 Hematopoietic cell lineage
  '
  dir="/home/bioinf/Desktop/newpaths"
  pv.out <- pathview(gene.data = genedata, pathway.id = pathwayid,
                   species = "human", out.suffix = c(n, pathwayid), kegg.dir=dir,
                   kegg.native = T)
}

for(pathwayid in list('05221', '05200', '05202', '04210', '04110', '04640')) {
  path_fun(et_1, 1, pathwayid)
  path_fun(et_12, 12, pathwayid)
  path_fun(et_6, 6, pathwayid)
}

# Example
data(gse16873.d)
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "05221",
                     species = "hsa", out.suffix = "gse16873")

# Make on for up and another for down

library(ggvenn)
x <- list(
  DE_1hr = row.names(et_1_toptags),
  DE_12hr = row.names(et_12_toptags),
  DE_6hr = row.names(et_6_toptags)
)


ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 3
)


