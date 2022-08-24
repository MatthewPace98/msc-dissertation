#######################################################################
# Combining DE genes from the treated samples into a single dataframe #
#######################################################################
# Assumes variables from edgeR.R are stored in memory

library(ggplot2)

# Extract counts for up and downregulated values
up_1hr = nrow(et_1_toptags[which(et_1_toptags$logFC > 0 ),])
down_1hr = nrow(et_1_toptags[which(et_1_toptags$logFC < 0 ),])
up_6hr = nrow(et_1_toptags[which(et_6_toptags$logFC > 0 ),])
down_6hr = nrow(et_1_toptags[which(et_6_toptags$logFC < 0 ),])
up_12hr = nrow(et_1_toptags[which(et_12_toptags$logFC > 0 ),])
down_12hr = nrow(et_1_toptags[which(et_12_toptags$logFC < 0 ),])


# Set variables for chart
time <- c('1hr', '12hr', '6hr', '1hr', '12hr', '6hr')
regulation <- c(rep('up', 3), rep('down', 3))
counts <- c(up_1hr,up_12hr,up_6hr,down_1hr,down_12hr,down_6hr)
DEG_counts <- data.frame(time, regulation, counts)
DEG_counts

colour = colorRampPalette(c("#0073C2FF", "#EFC000FF"), space = "rgb")(256)

# Open the drawing device.
png(file="/home/bioinf/Desktop/RNAseq/edgeR/DEG_counts_barchart.png")

# Draw the barchart
ggplot(DEG_counts, aes(fill=regulation, y=counts, x=time)) + 
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  scale_fill_manual( "", values = c("up" = "#0073C2FF", "down" = "#EFC000FF")) +
  labs(x = NULL, y = "Differentially Expressed Genes") +
  theme(text = element_text(size = 25)) + # increase text size
  scale_x_discrete(limits = c("1hr", "6hr", "12hr")) # Reorders x axis

# Turn off the device.
dev.off()

