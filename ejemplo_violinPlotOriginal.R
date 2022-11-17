load("~/Escritorio/SIC/TGFb2h_kmeans.RData")
load("~/Escritorio/SIC/TGFb2h_peaks_index.RData")
library(ggplot2)
library(ggpubr)
library(grid)

## Create index data 
cluster <- K_ATAC[[10]]$cluster

height <- c()
area <- c()
width.height.2 <- c()
local.peaks <- c()
M.height <- c()
clusterNum <- c()

for(i in 1:10){
  clust <- cluster==i
  height <- c(height, index.std_ATAC$height[clust])
  area <- c(area, index.std_ATAC$area[clust])
  width.height.2 <- c(width.height.2,index.std_ATAC$width.height.2[clust])
  local.peaks <- c(local.peaks,index.std_ATAC$local.peaks[clust])
  M.height <- c(M.height, index.std_ATAC$M.height[clust])
  clusterNum <- c(clusterNum, paste0("cluster_", rep(i, sum(clust))))
}

factor_clust = as.factor(clusterNum) ## Factor the number of clusters 
levels(factor_clust) = c(paste0(rep(1:10)))

featuresDF_ATAC<- data.frame("cluster"= factor_clust, height, area, width.height.2, local.peaks, M.height)

# Create the violin plots 
p1 <- ggplot(data=featuresDF_ATAC, aes(cluster, height)) + 
  geom_violin() + stat_summary(fun=mean, geom="point", size=2, color="red") + theme(axis.title.x=element_text(size = 15,hjust = 0.5),axis.title.y=element_blank(),axis.text.y = element_text(size = 15),axis.text.x = element_text(size = 15),plot.title = element_text(size=15,hjust = 0.5)) + ggtitle("Heigth")

p2 <-ggplot(data=featuresDF_ATAC, aes(cluster, area)) + 
  geom_violin() + stat_summary(fun=mean, geom="point", size=2, color="red") + theme(axis.title.x=element_text(size = 15,hjust = 0.5),axis.title.y=element_blank(),axis.text.y = element_text(size = 15),axis.text.x = element_text(size = 15),plot.title = element_text(size=15,hjust = 0.5)) + ggtitle("Area")

p3 <- ggplot(data=featuresDF_ATAC, aes(cluster, width.height.2)) + 
  geom_violin() + stat_summary(fun=mean, geom="point", size=2, color="red") + theme(axis.title.x=element_text(size = 15,hjust = 0.5),axis.title.y=element_blank(),axis.text.y = element_text(size = 15),axis.text.x = element_text(size = 15),plot.title = element_text(size=15,hjust = 0.5)) + ggtitle("Width at half maximun height")

p4 <-ggplot(data=featuresDF_ATAC, aes(cluster, local.peaks)) + 
  geom_violin() + stat_summary(fun=mean, geom="point", size=2, color="red") + theme(axis.title.x=element_text(size = 15,hjust = 0.5),axis.title.y=element_blank(),axis.text.y = element_text(size = 15),axis.text.x = element_text(size = 15),plot.title = element_text(size=15,hjust = 0.5)) + ggtitle("Local Peaks")

p5 <- ggplot(data=featuresDF_ATAC, aes(cluster, M.height)) + 
  geom_violin() + stat_summary(fun=mean, geom="point", size=2, color="red") + theme(axis.title.x=element_text(size = 15,hjust = 0.5),axis.title.y=element_blank(),axis.text.y = element_text(size = 15),axis.text.x = element_text(size = 15),plot.title = element_text(size =15,hjust = 0.5)) + ggtitle("M-Height")

# Save plots

jpeg("violin.jpeg",3000,500,units = "px")
ggarrange(p1,p2,p3,p4,p5,nrow = 1) 
dev.off()

#number of regions of each cluster
n_regions = c("Number of Cluster Regions",paste0("cluster_",rep(1:10)," = ",K[[10]][["size"]])) #save number of regions of each cluster
