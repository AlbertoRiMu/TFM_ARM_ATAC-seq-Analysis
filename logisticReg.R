setwd("Escritorio/SIC/results")

### Load indices data

load("TGFb2h_peaks.RData")
load("TGFb2hpeaks_index_std_ATAC.RData")


index <- data.frame("name" = peaks_ATAC$name, index.std_ATAC)


### Load region information (ENH OR PROM)
# The file we load is the regions that overlap with any TSS, information obtained from UCSC_mm9_knownGene
TSS <- read.table("TSS/ATAC_TX_TSS.bed")[,4]
prom <- c()
for(i in 1:length(peaks_ATAC$name)){
  if(sum(index$name[i]==TSS) == 1)
    prom[i] <- "prom"
  else
    prom[i] <- "enh"
}

index$promoter <- as.factor(prom)

### Variable distribution (INDICES)
# 1. Boxplot

jpeg("Index_by_prom.jpeg")
par(mfrow=c(2,3))
boxplot(height~promoter, index, outline= F, main = "Height",cex.main = 2,cex.axis = 1.4, cex.lab = 1.5,cex.names = 2)
boxplot(area~promoter, index, outline=F, main = "Area",cex.main = 2,cex.axis = 1.4, cex.lab = 1.5,cex.names = 2)
boxplot(width.height.2~promoter, index, outline=F, main = "Width.Height.2",cex.main = 2,cex.axis = 1.4, cex.lab = 1.5,cex.names = 2)
boxplot(local.peaks~promoter, index, outline=F, main= "local.peaks",cex.main = 2,cex.axis = 1.4, cex.lab = 1.5,cex.names = 2)
boxplot(M.height~promoter, index, outline=F, main="M.height",cex.main = 2,cex.axis = 1.4, cex.lab = 1.5,cex.names = 2)
dev.off()

# 2. collinearity

# pdf("logisticRegr/Index_correlation.pdf")
pairs(index[,c(2:6)])
# dev.off()

## Logistic regresion model
# Generating training and testing data set

set.seed(2022)

idx <- sample(nrow(index), nrow(index) * 0.25)
data.train <- index[idx, ] # Train the model
data.test <- index[-idx, ] # Test the model

# Generate model using glm and all the 5 variables

model <- glm(promoter~ height+ area  + width.height.2 + local.peaks +  M.height, 
             family = "binomial", data = data.train)
summary(model)

## Check for unneeded variables 
step(model, test="LRT")

### Testing the model. 
library(vcd)

predictions <- predict(model, newdata = data.test, type = "response")
dataPred <- data.frame(data.test$name, data.test$promoter, predictions)



predicciones <- ifelse(test = predictions > 0.3, yes = 1, no = 0) # Apply a 0.3 threshold to determine enh or prom
matriz_confusion <- table(ifelse(data.test$promoter =="prom", yes = 1, no = 0), 
                          predicciones, dnn = c("observaciones", "predicciones"))

# Confusion matrix
mosaic(matriz_confusion, shade = T, colorize = T, gp = gpar(fill = matrix(c("green3", "red2", "red2", "green3"), 2, 2)))
dataPredProm <- dataPred[dataPred$data.test.promoter=="prom",]

## PROMOTERS
matriz_confusion[2,2]/rowSums(matriz_confusion)[2]*100


## ENHANCERS
matriz_confusion[1,1]/rowSums(matriz_confusion)[1]*100

## Plot predictions vs true values
data.test2 <- data.frame(data.test, "odds" = predictions)
plot(odds~height, data.test2)
plot(odds~as.factor(local.peaks), data.test2)
plot(odds~area, data.test2)
par(mfrow= c(1,1))

## ROC Curve
library(pROC)
test_roc = roc(ifelse(data.test$promoter =="prom", yes = 1, no = 0) ~ predictions, plot = TRUE,
               print.auc = TRUE)
test_roc$auc
  
## Number of coincidences. Promoter regions in each cluster

prom_regions = index[index$promoter=="prom",]
prom_regions = prom_regions[,1]

intersect_1 = length(intersect(prom_regions,peaks1_ATAC))
intersect_2 = length(intersect(prom_regions,peaks2_ATAC))
intersect_3 = length(intersect(prom_regions,peaks3_ATAC))
intersect_4 = length(intersect(prom_regions,peaks4_ATAC))
intersect_5 = length(intersect(prom_regions,peaks5_ATAC))
intersect_6 = length(intersect(prom_regions,peaks6_ATAC))
intersect_7 = length(intersect(prom_regions,peaks7_ATAC))
intersect_8 = length(intersect(prom_regions,peaks8_ATAC))
intersect_9 = length(intersect(prom_regions,peaks9_ATAC))
intersect_10 = length(intersect(prom_regions,peaks10_ATAC))

vector_score = c(intersect_1,intersect_2,intersect_3,intersect_4,intersect_5,intersect_6,intersect_7,intersect_8,intersect_9,intersect_10)
vector_score_porcentaje_promoters = vector_score/length(index$name)*100
vector_names = c(1:10)

# Percentaje of each cluster over the prom regions
jpeg("% of each cluster over the promoter regions ATAC.jpeg")
barplot(vector_score_porcentaje_promoters,names.arg = vector_names,col = "royal blue",main = "% of clusters over promoter regions",xlab = "Clusters",ylab = "% of promoter regions",cex.main = 1.5,cex.axis = 1.2, cex.lab = 1.4,cex.names = 1.4)
dev.off()

# Percentaje of promoter regions on each cluster
vector_score_porcentaje_clusters = (vector_score/K_ATAC[[10]][["size"]])*100
jpeg("% of promoter regions per cluster ATAC.jpeg")
barplot(vector_score_porcentaje_clusters,names.arg = vector_names,col = "coral3",main = "% of promoters of each cluster",xlab = "Clusters",ylab = "% of promoter regions",cex.main = 1.5,cex.axis = 1.2, cex.lab = 1.4,cex.names = 1.4)
dev.off()


