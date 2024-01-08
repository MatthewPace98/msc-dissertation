#######################################################################
# Combining DE genes from the treated samples into a single dataframe #
#######################################################################
# Assumes variables from edgeR.R are stored in memory

library(org.Hs.eg.db)
library(AnnotationDbi) 
library(gplots)

# Takes the top 10 expressed genes of each time point (sorted by p value)
et_1_topten <- et_1_toptags
et_6_topten <- et_6_toptags
et_12_topten <- et_12_toptags

# Annotates genes with gene symbol
annotation_fun <- function(et){
  Symbol <- mapIds(org.Hs.eg.db, keys=rownames(et), keytype="ENSEMBL",
                   column="SYMBOL")
  data.frame(Symbol=Symbol)
}

# Renames rows to gene symbol names
et_1_topten$gene <- annotation_fun(et_1_topten)
row.names(et_1_topten) <- et_1_topten$gene$Symbol
et_12_topten$gene <- annotation_fun(et_12_topten)
row.names(et_12_topten) <- et_12_topten$gene$Symbol
et_6_topten$gene <- annotation_fun(et_6_topten)
row.names(et_6_topten) <- et_6_topten$gene$Symbol

# Extracts logFC columns of filtered datasets
temp_1 <- et_1_topten[c("logFC")]
colnames(temp_1)[colnames(temp_1) == 'logFC'] <- '1hr'

temp_6 <- et_6_topten[c("logFC")]
colnames(temp_6)[colnames(temp_6) == 'logFC'] <- '6hr'

temp_12 <- et_12_topten[c("logFC")]
colnames(temp_12)[colnames(temp_12) == 'logFC'] <- '12hr'

DE_combined <- merge(temp_1,temp_6,all=TRUE,by=0) # Merges dataframes
rownames(DE_combined)=DE_combined$Row.names
DE_combined <- DE_combined[-c(1)]

DE_combined <- merge(DE_combined, temp_12,all=TRUE,by=0) # Merges dataframes
rownames(DE_combined)=DE_combined$Row.names
DE_combined <- DE_combined[-c(1)]

DE_combined[is.na(DE_combined)] <- 0 # Changes NAs to 0

head(DE_combined)

# Check housekeeping genes
Symbol <- mapIds(org.Hs.eg.db, keys=Housekeeping_Genes$Ensembl, keytype="ENSEMBLTRANS",
                column="ENSEMBL")
Housekeeping_Genes$Ensembl.gene <- Symbol
head(Housekeeping_Genes$Ensembl.gene)
head(DE_combined)
intersect(rownames(DE_combined), Housekeeping_Genes$Ensembl.gene)

### Heatmap ### 
vals = as.matrix(DE_combined)

# Adds a little noise to each element to avoid the
# clustering function fail on zero variance datalines.
vals = jitter(vals, factor = 1, amount = 0.00001)

colour = colorRampPalette(c("#0073C2FF", "black", "#EFC000FF"), space = "rgb")(256)
png(file="/home/bioinf/Desktop/RNAseq/edgeR/heatmap_some.png")
heatmap.2(x=vals, 
          Colv=FALSE, 
          dendrogram="row",
          col=colour,
          trace="none",
          labRow=rownames(vals),
          main="Log2 Fold Changes of Top DEGs",
          ylab="")
dev.off()

