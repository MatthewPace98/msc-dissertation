library(org.Hs.eg.db)
library(AnnotationDbi) 
library(gplots)
library(gage)

# Assumes variables from edgeR.R are stored in memory

# Replaces Ensembl gene ID rownames in dataframe with Entrez IDs
entrez_fun <- function(et){
  logFC <- et[1] 
  logFC$entrez <- mapIds(org.Hs.eg.db, keys=rownames(et), 
                         keytype="ENSEMBL", column="ENTREZID")
  logFC <-na.omit(logFC)
  logFC <- logFC[duplicated(logFC[2]) == FALSE] 
  
  rownames(logFC) <- logFC$entrez
  logFC <- logFC[-c(2)]
  logFC
}
entrez_1 <- entrez_fun(et_1_toptags)
entrez_6 <- entrez_fun(et_6_toptags)
entrez_12 <- entrez_fun(et_12_toptags)

# Using GAGE to produce matrices of enriched gene sets and their respective statistics
gage_kegg_1 <- gage(entrez_1, gsets = kegg.gs, ref = NULL, samp = NULL)
gage_go_1 <- gage(entrez_1, gsets = go.gs, ref = NULL, samp = NULL)
gage_kegg_6 <- gage(entrez_6, gsets = kegg.gs, ref = NULL, samp = NULL)
gage_go_6 <- gage(entrez_6, gsets = go.gs, ref = NULL, samp = NULL)
gage_kegg_12 <- gage(entrez_12, gsets = kegg.gs, ref = NULL, samp = NULL)
gage_go_12 <- gage(entrez_12, gsets = go.gs, ref = NULL, samp = NULL)

# Removing gene sets with q<0.5 and NAs
gage_kegg_1 <- sigGeneSet(gage_kegg_1, cutoff = 0.5)
gage_go_1 <- sigGeneSet(gage_go_1, cutoff = 0.5)
gage_kegg_6 <- sigGeneSet(gage_kegg_6, cutoff = 0.5)
gage_go_6 <- sigGeneSet(gage_go_6, cutoff = 0.5)
gage_kegg_12 <- sigGeneSet(gage_kegg_12, cutoff = 0.5)
gage_go_12 <- sigGeneSet(gage_go_12, cutoff = 0.5)

# filters gene sets according to a q-value threshold and removes NAs
# sigGeneSet can create heatmaps but we opted out to keep our heatmaps standard
gage_kegg_1 <- sigGeneSet(gage_kegg_1, cutoff = 0.5)
gage_go_1 <- sigGeneSet(gage_go_1, cutoff = 0.5)
gage_kegg_6 <- sigGeneSet(gage_kegg_6, cutoff = 0.5)
gage_go_6 <- sigGeneSet(gage_go_6, cutoff = 0.5)
gage_kegg_12 <- sigGeneSet(gage_kegg_12, cutoff = 0.5)
gage_go_12 <- sigGeneSet(gage_go_12, cutoff = 0.5)

# Save as csv
path <-"/home/bioinf/Desktop/RNAseq/edgeR/GAGE/"
write.csv(rbind(gage_go_1$greater, gage_go_1$less), paste(path, "gage_go_1.csv"), row.names = TRUE)
write.csv(rbind(gage_kegg_1$greater, gage_kegg_1$less), paste(path, "gage_1_kegg.csv"), row.names = TRUE)
write.csv(rbind(gage_go_6$greater, gage_go_6$less), paste(path, "gage_go_6.csv"), row.names = TRUE)
write.csv(rbind(gage_kegg_6$greater, gage_kegg_6$less), paste(path, "gage_kegg_6.csv"), row.names = TRUE)
write.csv(rbind(gage_go_12$greater, gage_go_12$less), paste(path, "gage_go_12.csv"), row.names = TRUE)
write.csv(rbind(gage_kegg_12$greater, gage_kegg_12$less), paste(path, "gage_kegg_12.csv"), row.names = TRUE)

### Heatmap ### 
vals = as.matrix(gage_go_filt$stats[, 2:4])

# Adds a little noise to each element to avoid the
# clustering function fail on zero variance datalines.
vals = jitter(vals, factor = 1, amount = 0.00001)

colour = colorRampPalette(c("#0073C2FF", "black", "#EFC000FF"), space = "rgb")(256)
png(file="/home/bioinf/Desktop/RNAseq/edgeR/heatmap_gage_go.png", width = 700)
heatmap.2(x=vals, 
          Colv=FALSE, 
          dendrogram=NULL,
          col=colour,
          trace="none",
          labRow=rownames(vals),
          margins=c(5,35),
          main="Enriched GO Terms",
          ylab="")
dev.off()

### KEGG ###
vals = as.matrix(gage_kegg_filt$stats[, 2:4])

# Adds a little noise to each element to avoid the
# clustering function fail on zero variance datalines.
vals = jitter(vals, factor = 1, amount = 0.00001)

png(file="/home/bioinf/Desktop/RNAseq/edgeR/heatmap_gage_kegg.png", width = 570)
heatmap.2(x=vals, 
          Colv=FALSE, 
          dendrogram=NULL,
          col=colour,
          trace="none",
          labRow=rownames(vals),
          margins=c(5,21),
          main="Enriched KEGG pathways",
          ylab="")
dev.off()
