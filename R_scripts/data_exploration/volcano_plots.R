# Volcano plots
# Assumes variables from edgeR.R are stored in memory

draw <- function(et, et_toptags, time){
  # Open the drawing device.
  png(file=c("/home/bioinf/Desktop/RNAseq/edgeR/volcano/Volcano_6hr.png"), pointsize = 18)
  plot(et$table$logFC,-10*log10(et$table$PValue), 
       xlab="logFC", ylab="-10*log(p-value)", xlim=c(-10,10), ylim=c(0,1000))
  points(et_toptags$logFC, -10*log10(et_toptags$PValue), col="red", xlim=c(-10,10), ylim=c(0,1000))
  # Turn off the device.
  dev.off()
}

# Draw the Volcano plots
draw(et_1, et_1_toptags, "1hr")
draw(et_6, et_6_toptags, "6hr")
draw(et_12, et_12_toptags, "12hr")

