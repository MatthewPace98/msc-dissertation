# Venn diagram

library(ggvenn)

# To select for up and down regulated rows only, use something like:
#   '1hr' = row.names(et_1_toptags[which(et_1_toptags$logFC > 0 ),]),

x <- list(
  '1hr' = row.names(et_1_toptags),
  '12hr' = row.names(et_12_toptags),
  '6hr' = row.names(et_6_toptags)
)


png(file="/home/bioinf/Desktop/RNAseq/edgeR/venn.png")
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, 
  set_name_size = 6
)
dev.off()

