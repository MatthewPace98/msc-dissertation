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


# https://www.youtube.com/watch?v=8qvGNTVz3Ik
# input 'keys' i.e. genes.
# Specify the 'keytype', so ifthey are Entrez IDS, Ensembl etc
# Specify the output, so what you would like to convert them to
Symbol <- mapIds(org.Hs.eg.db, keys=rownames(y), keytype="ENSEMBL",
                 column=c("ENTREZID"))
y$genes <- data.frame(Symbol=Symbol)
head(y$genes)


# Filters genes with low counts 
# Threshold of 5-10 is normally taken
keep <- filterByExpr(y, group=group)
table(keep)
y <- y[keep,,keep.lib.sizes=FALSE]

# Normalisation by TMM
# y <- calcNormFactors(y)
design <- model.matrix(~group)

# Produces a plot in which distances between samples
# correspond to leading biological coefficient of variation (BCV) between those samples
#png(file="/home/bioinf/Desktop/RNAseq/edgeR/MDS.png")
#plotMDS(y)
#dev.off()

# Dispersion cannot be calculated so we take the BCV as 0.1 since we have no replicates
bcv <- 0.1
et <- exactTest(y, dispersion=bcv^2)

#Display top differentially expressed tags
# logFC, the log-abundance ratio, i.e. fold change, for each tag in the two groups being compared
# logCPM, the log-average concentration/abundance for each tag in the two groups being compared
# PValue, exact p-value for differential expression using the NB model
# FDR, the p-value adjusted for multiple testing as found using p.adjust using the method specified
classic_toptags <- topTags(et, n=20, adjust.method="BH", sort.by="logFC")$table

# To perform quasi-likelihood F-tests:
#QLfit <- glmQLFit(y,design, dispersion=bcv^2)
#qlf <- glmQLFTest(fit,coef=4,  dispersion=bcv^2)
#qlf_toptags <- topTags(qlf, n=20, adjust.method="BH", sort.by="logFC")$table

# To perform likelihood ratio tests:
fit <- glmFit(y,design, dispersion=bcv^2)
lrt <- glmLRT(fit,coef=4)
lrt_toptags <- topTags(lrt, n=20, adjust.method="BH", sort.by="logFC")$table

# Plot log-fold change against log-counts per million, with DE genes highlighted
#png(file="/home/bioinf/Desktop/RNAseq/edgeR/MD_4.png")
#plotMD(et)
#dev.off()


# A commonly used approach is to conduct DE tests, apply a fold-change cut-off 
# and then rank all the genes above that fold-change threshold by p-value
# The following is a rigorous statistical test for thresholded hypotheses.
# it tests whether the log2-fold-change is greater than lfc in absolute value.
#fit <- glmQLFit(y, design)
#tr <- glmTreat(fit, coef=ncol(fit), lfc=2)
#tr_toptags <- topTags(tr)$table


# KEGG and GO
geneid <-et$genes$Symbol

head(geneid)
go <- goana(et, species="Hs", geneid=geneid)
topGO(go, sort="up")
keg <- kegga(et, species="Hs", geneid=geneid)
topKEGG(keg, sort="up")

library(pathview)
# https://pathview.r-forge.r-project.org/pathview.pdf
# This package can be divided into four functional modules: the Downloader, 
# Parser, Mapper and Viewer
# gene.data can hold one sample or multiple sampleswith gene IDs as row names
# pathway.id the KEGG pathway ID usually 5 digits
# species is either the kegg code or name of target species. 'ko' is used for KEGG ortholog pathway
# kegg.native whether to render as png (TRUE) or graphiz layout (FASLE)
?pathview
data(gse16873.d)
data(demo.paths)
i <- 1
genedata <- as.matrix(et$table[, 1])
rownames(genedata) <- et$genes[, 1]
head(genedata)
pv.out <- pathview(gene.data = genedata, pathway.id = demo.paths$sel.paths[i],
                    species = "human", out.suffix = "gse16873_mine", kegg.native = TRUE)


# Heatmap
library(gplots)

# The gene names in the first column.
gene = row.names(y$counts)

# Load the data from the first column on.
vals = as.matrix(y$counts[, 1 : ncol(y$counts)])

# Adds a little noise to each element to avoid the
# clustering function fail on zero variance datalines.
vals = jitter(vals, factor = 1, amount = 0.00001)

# Each row is normalized to a z-score
zscore = NULL
for (i in 1 : nrow(vals)) {
  row = vals[i,]
  zrow = (row - mean(row)) / sd(row)
  zscore = rbind(zscore, zrow)
}

# Add back gene names as row names.
row.names(zscore) = gene

# Turn it into a matrix for heatmap2.
zscore = as.matrix(zscore)

# Open the drawing device.
pdf(file="/home/bioinf/Desktop/RNAseq/edgeR/heatmap.pdf")

# Set the color scheme.
colors = colorRampPalette(c("green", "black", "red"), space = "rgb")(256)

# Draw the heatmap.
# I think red is upregulated and green is down but need to check
heatmap.2(zscore, col = colors, density.info = "none", trace = "none", margins = c(7, 7), lhei = c(1, 5))

# Turn off the device.
dev.off()
