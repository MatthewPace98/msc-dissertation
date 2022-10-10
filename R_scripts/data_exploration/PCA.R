newdata <- subset(y$counts)
Group <- group
# Scree plot
pca<-prcomp(newdata,center=TRUE) 
plot(pca, type = "l") 
summary(pca) 
# 2 principal components explain 99.69% of the variability of the data 

# Biplot
library(devtools)
install_github("vqv/ggbiplot") ## downloading from github
library(ggbiplot)
biplot <- ggbiplot(pca, obs.scale = 1, var.scale = 1)  
biplot

autoplot(prcomp(newdata,center = TRUE), data = y$counts)

gbiplot <- ggbiplot(pca, obs.scale = 1, var.scale = 1, 
              groups = Group, ellipse = TRUE, 
              circle = TRUE)
gbiplot <- gbiplot + scale_color_discrete(name = '')
gbiplot <- gbiplot + theme(legend.direction = 'horizontal', 
               legend.position = 'bottom', aspect.ratio = 1)
gbiplot
ggsave(plot = gbiplot, width = 30, height = 30, dpi = 300, filename = "PCA.pdf")
