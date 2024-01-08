if (!requireNamespace("BiocManager", quietly=TRUE))
  + install.packages("BiocManager")
BiocManager::install(c("edgeR", "AnnotationDbi", "org.Hs.eg.db", "pathview", "gage", "gageData"))

library(edgeR)
library(pathview)
library(org.Hs.eg.db)
library(AnnotationDbi) 

files <- list("1hour.genes.results", "12hour.genes.results", "6hour.genes.results", "control.genes.results")
path <- "/home/bioinf/Desktop/RNAseq/RSEM_out/genes/"

# Creating a DGEList object from the estimated counts and gene identifier
labels = c('1hr', '12hr', '6hr', '0hr')
group <- factor(labels)
y <- readDGE(files, path=path, columns=c(1,5), group=group, labels=labels)

# Filter genes with low counts (<10)
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

# Differential Expression Analysis
# Dispersion cannot be calculated because we have no replicate so we assume the BCV is 0.1 
bcv <- 0.1
et_12 <- exactTest(y, pair=c(1, 2), dispersion=bcv^2)
et_1 <- exactTest(y, pair=c(1, 3), dispersion=bcv^2)
et_6 <- exactTest(y, pair=c(1, 4), dispersion=bcv^2)

#Display top differentially expressed tags
FDR_thresh <- 0.05
et_1_toptags <- topTags(et_1, n=nrow(et_1$table), adjust.method="BH", sort.by="PValue", p.value=FDR_thresh)$table
et_12_toptags <- topTags(et_12, n=nrow(et_12$table), adjust.method="BH", sort.by="PValue", p.value=FDR_thresh)$table
et_6_toptags <- topTags(et_6, n=nrow(et_6$table), adjust.method="BH", sort.by="PValue", p.value=FDR_thresh)$table

write.csv(y$counts,"/home/bioinf/Desktop/RNAseq/Top DE genes/filtered_counts.csv", row.names = TRUE)
write.csv(et_1_toptags,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/ettop1hr.csv", row.names = TRUE)
write.csv(et_12_toptags,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/ettop12hr.csv", row.names = TRUE)
write.csv(et_6_toptags,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/ettop6hr.csv", row.names = TRUE)

# Pathway analysis
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

# Look at pathways relevant to AML
pathway_ids <- list('05221', '05200', '05202', '04210', '04110', '04640', '04630', '04151', '04064', '04140', '04630')
for(pathwayid in pathway_ids) {
  path_fun(et_1_toptags, "1hr", pathwayid)
  path_fun(et_12_toptags, "12hr", pathwayid)
  path_fun(et_6_toptags, "6hr", pathwayid)
}

FC_thresh <- 1.5
et_1_toptags <- et_1_toptags[abs(et_1_toptags$logFC) > FC_thresh, ]
et_12_toptags <- et_12_toptags[abs(et_12_toptags$logFC) > FC_thresh, ]
et_6_toptags <- et_6_toptags[abs(et_6_toptags$logFC) > FC_thresh, ]

# Gene Set Enrichment Analysis
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


write.csv(et_1_toptags,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/ettop1hr.csv", row.names = TRUE)
write.csv(et_12_toptags,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/ettop12hr.csv", row.names = TRUE)
write.csv(et_6_toptags,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/ettop6hr.csv", row.names = TRUE)

write.csv(go_1,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/go_1hr.csv", row.names = TRUE)
write.csv(keg_1,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/kegg_1hr.csv", row.names = TRUE)

write.csv(go_12,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/go_12hr.csv", row.names = TRUE)
write.csv(keg_12,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/kegg_12hr.csv", row.names = TRUE)

write.csv(go_6,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/go_6hr.csv", row.names = TRUE)
write.csv(keg_6,"/home/bioinf/Desktop/RNAseq/edgeR/Top DE genes/kegg_6hr.csv", row.names = TRUE)
