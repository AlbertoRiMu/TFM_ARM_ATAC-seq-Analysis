 setwd("Escritorio/SIC/Matrix") 

# Code for getting matrix 
 # for i in TGFb2h_peak_class10  TGFb2h_peak_class1  TGFb2h_peak_class2  TGFb2h_peak_class3  TGFb2h_peak_class4  TGFb2h_peak_class5  TGFb2h_peak_class6  TGFb2h_peak_class7  TGFb2h_peak_class8  TGFb2h_peak_class9; do computeMatrix reference-point -p 20  --referencePoint center -a 3000 -b 3000 -bs 10 -R ../results/$i.bed -S ../input/matrix__TGFb2h.merged.bw --outFileName matrix/matrix.npz --outFileNameMatrix matrix/${i}.txt; done

# Create file list with the names and load them.
  files <- list.files(".", pattern = ".txt")
  names <- unlist(lapply(files, function(file) {gsub(file, pattern = "TGFb2h_peak_", replacement = "")}))
  names <- unlist(lapply(names, function(file) {gsub(file, pattern = ".txt", replacement = "")}))
  
  for(i in 1:length(names)){
    assign(names[i], read.table(files[i], header = T, sep = "\t"))
    }
  
  for(i in 1:length(names)){
    assign(names[i], apply(get(names[i]),2, mean, na.rm=T))
  }

## class03 = rep(NA,length(class03)) for cluster 3 of open data, the few reads mess with the scale 
dataF <- data.frame("bins" = rep(c(1:600), 10),
                      "coverage" = unlist(lapply(names, get)),
                      "class" = rep(names, each = 600))
  

# Plot matrix 

library(ggplot2)

jpeg("Matrix_Open_profile.jpeg", width = 70, height = 10, res = 200, units="cm")
ggplot(data=dataF, aes(x=bins, y=coverage, color= class))+
    geom_line() + facet_grid(~class) + theme_bw() + scale_x_continuous(name="bins", limits=c(0, 600),breaks = c(20,300,580),labels = c("-3kB","Center","3kB")) + theme(axis.text.x = element_text(size = 20,angle = 45,hjust = 1),axis.text.y =element_text(size=15),axis.title = element_text(size = 20),strip.text.x = element_text(size = 20))
dev.off()





