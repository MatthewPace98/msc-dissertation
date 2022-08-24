library(org.Hs.eg.db)
library(AnnotationDbi) 

# Assumes variables from edgeR.R are stored in memory
library(gage)

# GAGE uses Entrez IDs
# Replaces Ensembl gene ID rownames in dataframe with Entrez IDs
entrez_fun <- function(et){
  # Input Ensembl ID matrix of DEGs, logFCs, etc.
  # Output Entrez ID matrix of DEGs and logFC
  logFC <- et[1] 
  logFC$entrez <- mapIds(org.Hs.eg.db, keys=rownames(et), 
                         keytype="ENSEMBL", column="ENTREZID") # create EntrezID column
  logFC <-na.omit(logFC)
  logFC <- logFC[duplicated(logFC[2]) == FALSE] 
  
  rownames(logFC) <- logFC$entrez # set EntrezID as rows
  logFC <- logFC[-c(2)] # remove EntrezID column
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
# The default cutoff is 0.1
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

###############
### Heatmap ### 
###############
library(gplots)

# Load the data 
vals = as.matrix(gage_go_filt$stats[, 2:4])

# Adds a little noise to each element to avoid the
# clustering function fail on zero variance datalines.
vals = jitter(vals, factor = 1, amount = 0.00001)

# Set the color scheme.
colour = colorRampPalette(c("#0073C2FF", "black", "#EFC000FF"), space = "rgb")(256)

# Open the drawing device.
png(file="/home/bioinf/Desktop/RNAseq/edgeR/heatmap_gage_go.png", width = 700)

# Draw the heatmap.
heatmap.2(x=vals, 
          Colv=FALSE, 
          dendrogram=NULL,
          col=colour,
          trace="none",
          labRow=rownames(vals),
          margins=c(5,35),
          main="Enriched GO Terms",
          ylab="")

# Turn off the device.
dev.off()

### KEGG ###
# Load the data 
vals = as.matrix(gage_kegg_filt$stats[, 2:4])

# Adds a little noise to each element to avoid the
# clustering function fail on zero variance datalines.
vals = jitter(vals, factor = 1, amount = 0.00001)

# Open the drawing device.
png(file="/home/bioinf/Desktop/RNAseq/edgeR/heatmap_gage_kegg.png", width = 570)

# Draw the heatmap.
# heatmap with the defaults parameters
heatmap.2(x=vals, 
          Colv=FALSE, 
          dendrogram=NULL,
          col=colour,
          trace="none",
          labRow=rownames(vals),
          margins=c(5,21),
          main="Enriched KEGG pathways",
          ylab="")

# Turn off the device.
dev.off()
