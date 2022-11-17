##### Heatmaps #####

# Load paths to data
path_ATAC = "results_ATAC/TGFb2h_peak_class"
path_Monoc= "results_monoc/TGFb2h_peak_class"
path_open = "results_open/TGFb2h_peak_class"
path_Monoc_open = "results_Monoc_Open/TGFb2h_peak_class"

# Load cluster data. Repeat from 1 to 10
clust1_ATAC= read.table(paste0(path_ATAC,1,".bed"),header = F, sep = "\t")
clust1_Monoc_open = read.table(paste0(path_Monoc_open,1,".bed"),header = F, sep = "\t")
clust1_monoc= read.table(paste0(path_Monoc,1,".bed"),header = F, sep = "\t")
clust1_open= read.table(paste0(path_open,1,".bed"),header = F, sep = "\t")

# Get the names of hte peaks. Repeat from 1 to 10
colnames(clust1) <- c("chr","start","end","peak","score","strand")
colnames(clust1) <- c("chr","start","end","peak","score","strand")
peaks1_ATAC<- clust1_ATAC$peak
peaks1_Monoc<- clust1_monoc$peak
peaks1_open <- clust1_open$peak
peaks10_Monoc_open <- clust1_Monoc_open$peak
 

## Counting of intersections. Example with Mononuc and Open
library(vecsets)

score_Monoc_open_1_1 = length(intersect(peaks1_open,peaks1_Monoc))
score_Monoc_open_1_2 = length(intersect(peaks1_open,peaks2_Monoc))
score_Monoc_open_1_3 = length(intersect(peaks1_open,peaks3_Monoc))
score_Monoc_open_1_4 = length(intersect(peaks1_open,peaks4_Monoc))
score_Monoc_open_1_5 = length(intersect(peaks1_open,peaks5_Monoc))
score_Monoc_open_1_6 = length(intersect(peaks1_open,peaks6_Monoc))
score_Monoc_open_1_7 = length(intersect(peaks1_open,peaks7_Monoc))
score_Monoc_open_1_8 = length(intersect(peaks1_open,peaks8_Monoc))
score_Monoc_open_1_9 = length(intersect(peaks1_open,peaks9_Monoc))
score_Monoc_open_1_10 = length(intersect(peaks1_open,peaks10_Monoc))
score_Monoc_open_2_1 = length(intersect(peaks2_open,peaks1_Monoc))
score_Monoc_open_2_2 = length(intersect(peaks2_open,peaks2_Monoc))
score_Monoc_open_2_3 = length(intersect(peaks2_open,peaks3_Monoc))
score_Monoc_open_2_4 = length(intersect(peaks2_open,peaks4_Monoc))
score_Monoc_open_2_5 = length(intersect(peaks2_open,peaks5_Monoc))
score_Monoc_open_2_6 = length(intersect(peaks2_open,peaks6_Monoc))
score_Monoc_open_2_7 = length(intersect(peaks2_open,peaks7_Monoc))
score_Monoc_open_2_8 = length(intersect(peaks2_open,peaks8_Monoc))
score_Monoc_open_2_9 = length(intersect(peaks2_open,peaks9_Monoc))
score_Monoc_open_2_10 = length(intersect(peaks2_open,peaks10_Monoc))
score_Monoc_open_3_1 = length(intersect(peaks3_open,peaks1_Monoc))
score_Monoc_open_3_2 = length(intersect(peaks3_open,peaks2_Monoc))
score_Monoc_open_3_3 = length(intersect(peaks3_open,peaks3_Monoc))
score_Monoc_open_3_4 = length(intersect(peaks3_open,peaks4_Monoc))
score_Monoc_open_3_5 = length(intersect(peaks3_open,peaks5_Monoc))
score_Monoc_open_3_6 = length(intersect(peaks3_open,peaks6_Monoc))
score_Monoc_open_3_7 = length(intersect(peaks3_open,peaks7_Monoc))
score_Monoc_open_3_8 = length(intersect(peaks3_open,peaks8_Monoc))
score_Monoc_open_3_9 = length(intersect(peaks3_open,peaks9_Monoc))
score_Monoc_open_3_10 = length(intersect(peaks3_open,peaks10_Monoc))
score_Monoc_open_4_1 = length(intersect(peaks4_open,peaks1_Monoc))
score_Monoc_open_4_2 = length(intersect(peaks4_open,peaks2_Monoc))
score_Monoc_open_4_3 = length(intersect(peaks4_open,peaks3_Monoc))
score_Monoc_open_4_4 = length(intersect(peaks4_open,peaks4_Monoc))
score_Monoc_open_4_5 = length(intersect(peaks4_open,peaks5_Monoc))
score_Monoc_open_4_6 = length(intersect(peaks4_open,peaks6_Monoc))
score_Monoc_open_4_7 = length(intersect(peaks4_open,peaks7_Monoc))
score_Monoc_open_4_8 = length(intersect(peaks4_open,peaks8_Monoc))
score_Monoc_open_4_9 = length(intersect(peaks4_open,peaks9_Monoc))
score_Monoc_open_4_10 = length(intersect(peaks4_open,peaks10_Monoc))
score_Monoc_open_5_1 = length(intersect(peaks5_open,peaks1_Monoc))
score_Monoc_open_5_2 = length(intersect(peaks5_open,peaks2_Monoc))
score_Monoc_open_5_3 = length(intersect(peaks5_open,peaks3_Monoc))
score_Monoc_open_5_4 = length(intersect(peaks5_open,peaks4_Monoc))
score_Monoc_open_5_5 = length(intersect(peaks5_open,peaks5_Monoc))
score_Monoc_open_5_6 = length(intersect(peaks5_open,peaks6_Monoc))
score_Monoc_open_5_7 = length(intersect(peaks5_open,peaks7_Monoc))
score_Monoc_open_5_8 = length(intersect(peaks5_open,peaks8_Monoc))
score_Monoc_open_5_9 = length(intersect(peaks5_open,peaks9_Monoc))
score_Monoc_open_5_10 = length(intersect(peaks5_open,peaks10_Monoc))
score_Monoc_open_6_1 = length(intersect(peaks6_open,peaks1_Monoc))
score_Monoc_open_6_2 = length(intersect(peaks6_open,peaks2_Monoc))
score_Monoc_open_6_3 = length(intersect(peaks6_open,peaks3_Monoc))
score_Monoc_open_6_4 = length(intersect(peaks6_open,peaks4_Monoc))
score_Monoc_open_6_5 = length(intersect(peaks6_open,peaks5_Monoc))
score_Monoc_open_6_6 = length(intersect(peaks6_open,peaks6_Monoc))
score_Monoc_open_6_7 = length(intersect(peaks6_open,peaks7_Monoc))
score_Monoc_open_6_8 = length(intersect(peaks6_open,peaks8_Monoc))
score_Monoc_open_6_9 = length(intersect(peaks6_open,peaks9_Monoc))
score_Monoc_open_6_10 = length(intersect(peaks6_open,peaks10_Monoc))
score_Monoc_open_7_1 = length(intersect(peaks7_open,peaks1_Monoc))
score_Monoc_open_7_2 = length(intersect(peaks7_open,peaks2_Monoc))
score_Monoc_open_7_3 = length(intersect(peaks7_open,peaks3_Monoc))
score_Monoc_open_7_4 = length(intersect(peaks7_open,peaks4_Monoc))
score_Monoc_open_7_5 = length(intersect(peaks7_open,peaks5_Monoc))
score_Monoc_open_7_6 = length(intersect(peaks7_open,peaks6_Monoc))
score_Monoc_open_7_7 = length(intersect(peaks7_open,peaks7_Monoc))
score_Monoc_open_7_8 = length(intersect(peaks7_open,peaks8_Monoc))
score_Monoc_open_7_9 = length(intersect(peaks7_open,peaks9_Monoc))
score_Monoc_open_7_10 = length(intersect(peaks7_open,peaks10_Monoc))
score_Monoc_open_8_1 = length(intersect(peaks8_open,peaks1_Monoc))
score_Monoc_open_8_2 = length(intersect(peaks8_open,peaks2_Monoc))
score_Monoc_open_8_3 = length(intersect(peaks8_open,peaks3_Monoc))
score_Monoc_open_8_4 = length(intersect(peaks8_open,peaks4_Monoc))
score_Monoc_open_8_5 = length(intersect(peaks8_open,peaks5_Monoc))
score_Monoc_open_8_6 = length(intersect(peaks8_open,peaks6_Monoc))
score_Monoc_open_8_7 = length(intersect(peaks8_open,peaks7_Monoc))
score_Monoc_open_8_8 = length(intersect(peaks8_open,peaks8_Monoc))
score_Monoc_open_8_9 = length(intersect(peaks8_open,peaks9_Monoc))
score_Monoc_open_8_10 = length(intersect(peaks8_open,peaks10_Monoc))
score_Monoc_open_9_1 = length(intersect(peaks9_open,peaks1_Monoc))
score_Monoc_open_9_2 = length(intersect(peaks9_open,peaks2_Monoc))
score_Monoc_open_9_3 = length(intersect(peaks9_open,peaks3_Monoc))
score_Monoc_open_9_4 = length(intersect(peaks9_open,peaks4_Monoc))
score_Monoc_open_9_5 = length(intersect(peaks9_open,peaks5_Monoc))
score_Monoc_open_9_6 = length(intersect(peaks9_open,peaks6_Monoc))
score_Monoc_open_9_7 = length(intersect(peaks9_open,peaks7_Monoc))
score_Monoc_open_9_8 = length(intersect(peaks9_open,peaks8_Monoc))
score_Monoc_open_9_9 = length(intersect(peaks9_open,peaks9_Monoc))
score_Monoc_open_9_10 = length(intersect(peaks9_open,peaks10_Monoc))
score_Monoc_open_10_1 = length(intersect(peaks10_open,peaks1_Monoc))
score_Monoc_open_10_2 = length(intersect(peaks10_open,peaks2_Monoc))
score_Monoc_open_10_3 = length(intersect(peaks10_open,peaks3_Monoc))
score_Monoc_open_10_4 = length(intersect(peaks10_open,peaks4_Monoc))
score_Monoc_open_10_5 = length(intersect(peaks10_open,peaks5_Monoc))
score_Monoc_open_10_6 = length(intersect(peaks10_open,peaks6_Monoc))
score_Monoc_open_10_7 = length(intersect(peaks10_open,peaks7_Monoc))
score_Monoc_open_10_8 = length(intersect(peaks10_open,peaks8_Monoc))
score_Monoc_open_10_9 = length(intersect(peaks10_open,peaks9_Monoc))
score_Monoc_open_10_10 = length(intersect(peaks10_open,peaks10_Monoc))

vector_score_open_Monoc= c(score_Monoc_open_1_1 ,  score_Monoc_open_1_2 , score_Monoc_open_1_3 , score_Monoc_open_1_4 , score_Monoc_open_1_5 ,  score_Monoc_open_1_6 , score_Monoc_open_1_7 , score_Monoc_open_1_8 , score_Monoc_open_1_9 , score_Monoc_open_1_10 , score_Monoc_open_2_1 , score_Monoc_open_2_2 , score_Monoc_open_2_3 , score_Monoc_open_2_4 , score_Monoc_open_2_5 , score_Monoc_open_2_6 , score_Monoc_open_2_7 , score_Monoc_open_2_8 , score_Monoc_open_2_9 , score_Monoc_open_2_10 , score_Monoc_open_3_1 , score_Monoc_open_3_2 , score_Monoc_open_3_3 , score_Monoc_open_3_4 , score_Monoc_open_3_5 , score_Monoc_open_3_6 , score_Monoc_open_3_7 , score_Monoc_open_3_8 ,  score_Monoc_open_3_9 , score_Monoc_open_3_10 , score_Monoc_open_4_1 ,  score_Monoc_open_4_2 , score_Monoc_open_4_3 , score_Monoc_open_4_4 , score_Monoc_open_4_5 , score_Monoc_open_4_6 , score_Monoc_open_4_7 , score_Monoc_open_4_8 , score_Monoc_open_4_9 , score_Monoc_open_4_10 , score_Monoc_open_5_1 , score_Monoc_open_5_2 , score_Monoc_open_5_3 , score_Monoc_open_5_4 , score_Monoc_open_5_5 , score_Monoc_open_5_6 , score_Monoc_open_5_7 , score_Monoc_open_5_8 , score_Monoc_open_5_9 , score_Monoc_open_5_10 , score_Monoc_open_6_1 , score_Monoc_open_6_2 , score_Monoc_open_6_3 , score_Monoc_open_6_4 , score_Monoc_open_6_5 , score_Monoc_open_6_6 , score_Monoc_open_6_7 , score_Monoc_open_6_8 , score_Monoc_open_6_9 , score_Monoc_open_6_10 , score_Monoc_open_7_1 , score_Monoc_open_7_2 , score_Monoc_open_7_3 , score_Monoc_open_7_4 , score_Monoc_open_7_5 , score_Monoc_open_7_6 , score_Monoc_open_7_7 , score_Monoc_open_7_8 ,  score_Monoc_open_7_9 , score_Monoc_open_7_10 ,  score_Monoc_open_8_1 ,  score_Monoc_open_8_2 ,  score_Monoc_open_8_3 ,  score_Monoc_open_8_4 ,  score_Monoc_open_8_5 ,  score_Monoc_open_8_6 , score_Monoc_open_8_7 ,  score_Monoc_open_8_8 ,  score_Monoc_open_8_9 ,  score_Monoc_open_8_10 ,  score_Monoc_open_9_1 ,  score_Monoc_open_9_2 ,  score_Monoc_open_9_3 , score_Monoc_open_9_4 ,  score_Monoc_open_9_5 ,  score_Monoc_open_9_6 ,  score_Monoc_open_9_7 , score_Monoc_open_9_8 ,  score_Monoc_open_9_9 , score_Monoc_open_9_10 , score_Monoc_open_10_1 ,  score_Monoc_open_10_2 ,  score_Monoc_open_10_3 , score_Monoc_open_10_4 ,  score_Monoc_open_10_5  ,  score_Monoc_open_10_6 ,  score_Monoc_open_10_7 ,  score_Monoc_open_10_8 ,  score_Monoc_open_10_9 ,  score_Monoc_open_10_10)


matrix_score= matrix(vector_score_ATAC_Monoc_open,ncol = 10, nrow = 10)

# Check if rows or columns sum total of monoc or open regions
sum(matrix_score[,1]) ==  K_mononuc[[10]][["size"]][1] ## False
sum(matrix_score[1,]) ==  K_mononuc[[10]][["size"]][1] ## True. Mononuc regions are the sum of rows

# Heatmap Column scaled
matrix2_score_ATAC_Monoc_open_percentaje = matrix(rep(NA),100,ncol = 10, nrow = 10)
for (i in 1:nrow(matrix_score)){
  for (j in 1:ncol(matrix_score)){
    matrix2_score_ATAC_Monoc_open_percentaje[i,j] = (matrix_score[i,j]/sum(matrix_score[,j]))*100
  }
  
}

colnames(matrix2_score_ATAC_Monoc_open_percentaje) <- c(paste0(seq(1:10)))
rownames(matrix2_score_ATAC_Monoc_open_percentaje) <- c(paste0(seq(1:10)))
out_heatmap = "Heatmaps/Heatmap_"


## Heatmap row scaled

matrix2_score_percentaje_R = matrix(rep(NA),100,ncol = 10, nrow = 10)
for (i in 1:nrow(matrix_score)){
  for (j in 1:ncol(matrix_score)){
    matrix2_score_percentaje_R[i,j] = (matrix_score[i,j]/sum(matrix_score[i,]))*100
  }
  
}

colnames(matrix2_score_percentaje_R) <- c(paste0(seq(1:10)))
rownames(matrix2_score_percentaje_R) <- c(paste0(seq(1:10)))


#Creating heatmaps 
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
vector = seq(0,100,length.out = 9)
col_fun = colorRamp2(c(0,25,50,75, 100),brewer.pal(5, "YlOrRd"))

# Col scaled
Heatmap(matrix2_score_ATAC_Monoc_open_percentaje, cluster_rows = F,cluster_columns = F, col = col_fun, column_names_gp = grid::gpar(fontsize = 30),
  row_names_gp = grid::gpar(fontsize = 30))


# Row scaled
Heatmap(matrix2_score_percentaje_R, cluster_rows = F,cluster_columns = F, col = col_fun, heatmap_legend_param = list(title = "% of Coincidence Mononuc scaled",legend_height = unit(6, "cm"),grid_width = unit(2, "cm"),title_gp=gpar(fontsize=12, fontface="bold"),labels_gp = gpar(fontsize = 14)))





, heatmap_legend_param = list(title = "% of Coincidence ATAC scaled",legend_height = unit(8, "cm"),grid_width = unit(2, "cm"),title_gp=gpar(fontsize=16, fontface="bold"),labels_gp = gpar(fontsize = 14))
