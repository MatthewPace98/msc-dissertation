# Get log2 counts
logcounts <- log2(y$counts + 1)

# Check distributions of samples using boxplots
png(file="/home/bioinf/Desktop/RNAseq/edgeR/boxplot_filtered.png")
boxplot(logcounts,
        xlab="",
        ylab="Log2(Counts)",
        las=2,
        main="Filtered Log2(Counts)")

# Add a blue line that corresponds to the median
abline(h=median(logcounts), col="blue")
dev.off()