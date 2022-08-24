if (!requireNamespace("BiocManager", quietly=TRUE))
  + install.packages("BiocManager")
BiocManager::install(c("edgeR", "AnnotationDbi", "org.Hs.eg.db", "pathview", "gage", "gageData"))

library(edgeR)

# https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# https://www.bioconductor.org/packages/devel/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html
# https://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf

files <- list("1hour.genes.results", "12hour.genes.results", "6hour.genes.results", "control.genes.results")
path <- "/home/bioinf/Desktop/RNAseq/RSEM_out/genes/"

# Creating a DGEList object from the estimated counts and gene identifier
labels = c('1hr', '12hr', '6hr', '0hr')
group <- factor(labels)
y <- readDGE(files, path=path, columns=c(1,5), group=group, labels=labels)

# Filters genes with low counts (<10)
keep <- filterByExpr(y, group=group)
table(keep)  # how many were filtered and how many remain
y <- y[keep,keep.lib.sizes=FALSE]

# Normalisation by TMM
y <- calcNormFactors(y)
design <- model.matrix(~group)


# Produces a PCA plot
png(file="/home/bioinf/Desktop/RNAseq/edgeR/PCA.png", pointsize = 16)
plotMDS(y, method='logFC', gene.selection='common')
dev.off()

# Dispersion cannot be calculated so we take the BCV as 0.1 since we have no replicates
bcv <- 0.1
et_12 <- exactTest(y, pair=c(1, 2), dispersion=bcv^2)
et_1 <- exactTest(y, pair=c(1, 3), dispersion=bcv^2)
et_6 <- exactTest(y, pair=c(1, 4), dispersion=bcv^2)

# Plot log-fold change against log-counts per million, with DE genes highlighted
# What is this doing? Maybe another volcano plot
#png(file="/home/bioinf/Desktop/RNAseq/edgeR/MD/MD_1.png")
#plotMD(et_1, column = 1)
#dev.off()
#png(file="/home/bioinf/Desktop/RNAseq/edgeR/MD/MD_12.png")
#plotMD(et_12)
#dev.off()
#png(file="/home/bioinf/Desktop/RNAseq/edgeR/MD/MD_6.png")
#plotMD(et_6)
#dev.off()

#Display top differentially expressed tags
# logFC, the log2 fold change, for each tag in the two groups being compared
# logCPM, the log-average concentration/abundance for each tag in the two groups being compared
# PValue, p-value for differential expression using the NB model
# FDR, the p-value adjusted for multiple testing as found using p.adjust using the method specified
# p.value is filtering by adjusted p (i.e. FDR), not raw p 
FDR_thresh <- 0.05  # filters rows with FDR less than this
et_1_toptags <- topTags(et_1, n=nrow(et_1$table), adjust.method="BH", sort.by="PValue", p.value=FDR_thresh)$table
et_12_toptags <- topTags(et_12, n=nrow(et_12$table), adjust.method="BH", sort.by="PValue", p.value=FDR_thresh)$table
et_6_toptags <- topTags(et_6, n=nrow(et_6$table), adjust.method="BH", sort.by="PValue", p.value=FDR_thresh)$table

write.csv(y$counts,"/home/bioinf/Desktop/RNAseq/Top DE genes/filtered_counts.csv", row.names = TRUE)
write.csv(et_1_toptags,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/ettop1hr.csv", row.names = TRUE)
write.csv(et_12_toptags,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/ettop12hr.csv", row.names = TRUE)
write.csv(et_6_toptags,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/ettop6hr.csv", row.names = TRUE)

library(pathview)
# https://pathview.r-forge.r-project.org/pathview.pdf
# This package can be divided into four functional modules: the Downloader, 
# Parser, Mapper and Viewer
# pathway.id is the KEGG pathway ID
# kegg.native whether to render as png (TRUE) or graphiz layout (FALSE)
# KEGG pathways: https://www.genome.jp/kegg/pathway.html
# No difference between Entrez and Ensembl annotations in erms of pathways
# Should this come before p value cutoffs? I dont think so
if (1 == 2){
  '
  05221 Acute myeloid leukemia
  05200 Pathways in cancer
  05202 Transcriptional misregulation in cancer
  04210 Apoptosis
  04110 Cell cycle
  04640 Hematopoietic cell lineage
  04630 JAK-STAT signaling pathway
  04151 PI3K-Akt signaling pathway
  04064 NF-kappa B pathway
  04140 Autophagy
  04630 JAK 
'
  
  path_fun <- function(et, sample, pathwayid){
    genedata <- as.matrix(et$logFC)  
    rownames(genedata) <- rownames(et)
    head(genedata)
    dir="/home/bioinf/Desktop/RNAseq/edgeR/pathview"
    
    download.kegg(pathway.id = pathwayid, species = "hsa", kegg.dir = dir,
                  file.type=c("xml", "png"))
    pv.out <- pathview(gene.data = genedata, pathway.id = pathwayid,
                       species = "human", out.suffix = c(sample, pathwayid), kegg.dir=dir,
                       gene.idtype = "ENSEMBL", kegg.native = T)
  }
  
  
  for(pathwayid in list('05221', '05200', '05202', '04210', '04110', '04640', '04630', '04151', '04064', '04140', '04630')) {
    path_fun(et_1_toptags, "1hr", pathwayid)
    path_fun(et_12_toptags, "12hr", pathwayid)
    path_fun(et_6_toptags, "6hr", pathwayid)
  }
}
FC_thresh <- 1.5  # filters rows with a logFC less than this
et_1_toptags <- et_1_toptags[abs(et_1_toptags$logFC) > FC_thresh, ]
et_12_toptags <- et_12_toptags[abs(et_12_toptags$logFC) > FC_thresh, ]
et_6_toptags <- et_6_toptags[abs(et_6_toptags$logFC) > FC_thresh, ]

# Replaces et with et_toptags
# et_1$table <- et_1_toptags[2:4]
# et_1$genes <- et_1_toptags[1]

# et_12$table <- et_12_toptags[2:4]
# et_12$genes <- et_12_toptags[1]

# et_6$table <- et_6_toptags[2:4]
# et_6$genes <- et_6_toptags[1]

#########################################
### Gene Set testing with KEGG and GO ###
#########################################

library(org.Hs.eg.db)
library(AnnotationDbi) 
# 'keys' = genes
# 'keytype' = type of ID (entrez, ensembl, etc)
# 'column' = output format ("SYMBOL", "GO", "ENTREZID")
annotation_fun <- function(et){
  Symbol <- mapIds(org.Hs.eg.db, keys=rownames(et), keytype="ENSEMBL",
                   column="ENTREZID")
  et$Symbol <- data.frame(Symbol=Symbol)
  et$Symbol
}

# Geneid needs to be Entrez IDs
geneid_12 <-annotation_fun(et_12_toptags)
head(geneid_12)
colMeans(is.na(geneid_12))  # Proportion of NAs

geneid_1 <-annotation_fun(et_1_toptags)
head(geneid_1)
colMeans(is.na(geneid_1))

geneid_6 <-annotation_fun(et_6_toptags)
head(geneid_6)
colMeans(is.na(geneid_6))

go_1 <- goana(geneid_1$Symbol, species="Hs", geneid=geneid_1)
keg_1 <- kegga(geneid_1$Symbol, species="Hs", geneid=geneid_1)

go_12 <- goana(geneid_12$Symbol, species="Hs", geneid=geneid_12)
keg_12 <- kegga(geneid_12$Symbol, species="Hs", geneid=geneid_12)

go_6 <- goana(geneid_6$Symbol, species="Hs", geneid=geneid_6)
keg_6 <- kegga(geneid_6$Symbol, species="Hs", geneid=geneid_6)



# write.csv(y$counts,"/home/bioinf/Desktop/RNAseq/Top DE genes/filtered_counts.csv", row.names = TRUE)
# write.csv(lrt_toptags,"/home/bioinf/Desktop/RNAseq/Top DE genes/lrttop.csv", row.names = TRUE)
write.csv(et_1_toptags,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/ettop1hr.csv", row.names = TRUE)
write.csv(et_12_toptags,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/ettop12hr.csv", row.names = TRUE)
write.csv(et_6_toptags,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/ettop6hr.csv", row.names = TRUE)

write.csv(go_1,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/go_1hr.csv", row.names = TRUE)
write.csv(keg_1,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/kegg_1hr.csv", row.names = TRUE)

write.csv(go_12,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/go_12hr.csv", row.names = TRUE)
write.csv(keg_12,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/kegg_12hr.csv", row.names = TRUE)

write.csv(go_6,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/go_6hr.csv", row.names = TRUE)
write.csv(keg_6,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/kegg_6hr.csv", row.names = TRUE)
