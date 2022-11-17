require(methods)
require(rtracklayer)
require(Rcpp)
setwd("~/Escritorio/SIC/")

###############################
##### PEAK IDENTIFICATION #####
###############################

importPeaks.bw2 <- function(peaksRange,count.bw,N="all",start=1){
  require(rtracklayer)
  peaks <- import(peaksRange)
  if(N!="all" & N<length(peaks)){
    end=min(length(peaks),start+N-1)
    peaks=peaks[start:end]
  }
  n_peaks=length(peaks)
  
  library(parallel) 
  core.number <- max(min(floor(n_peaks/10),detectCores()),1) # at least 10 peaks for each node
  
  importPeak=function(peak,count.bw){
    count=import(count.bw,format="bw",which=peak)
    scores=as.vector(Rle(count$score,width(count)))
    if(sum(width(count))!=(end(count[length(count)])-start(count[1])+1)){
      gaps=c(start(count[-1])-end(count[-length(count)])-1,0)
      index=cumsum(width(count))[gaps!=0]
      for(j in length(index):1)
        scores=c(scores[1:index[j]],rep(0,gaps[gaps!=0][j]),scores[(index[j]+1):length(scores)])
    }
    diff_start=start(peak)-start(count[1])
    if(diff_start<0){
      scores=c(rep(0,times=-diff_start),scores)
      diff_start=0
    }
    count=scores[(diff_start+1):(diff_start+width(peak))]
    count[is.na(count)]=0
    return(count)
  }
  n_peaks_group=min(floor(n_peaks/core.number),100)
  groups=c(rep.int(1,n_peaks-floor(n_peaks/n_peaks_group)*n_peaks_group),rep(seq_len(floor(n_peaks/n_peaks_group)),each=n_peaks_group))
  cl <- makeCluster(core.number)
  peaks$count=Reduce(c,parLapplyLB(cl,split(peaks,groups),
                                   function(peaks_group,count.bw,importPeak){
                                     require(rtracklayer)
                                     return(lapply(seq_along(peaks_group),function(i) importPeak(peaks_group[i],count.bw)))
                                   },count.bw=count.bw,importPeak=importPeak))
  stopCluster(cl)
  return(peaks)
}



##### Prepare files
dir.create("results/")

bw= "input/ATAC_TGFb2h_MonoNuc.bw" # this will be done with ATAC_TGFb2h.bw and ATAC_TGFb2h_open.bw too
bed= "input/atac_seq_peaks.bed"
out= "results/TGFb2h"


# LOAD FUNCTIONS FOR PEAK CHARACTERISTICS
load(paste0("~/opt/SIC-ChIP-master/SIC-ChIP_functions.RData"))

# COMPILADOR C++ code
Rcpp::sourceCpp("scripts/picco.cpp")

#######################################################
#### import coverage function for each peak in bed ####
#######################################################
peaks <- importPeaks.bw2(bed,bw)

plot(c(1:length(unlist(peaks$count[1000]))), unlist(peaks$count[1000]), type = "l")

removePeaks <- c()
for(i in 1:length(peaks$name)){
  removePeaks[i] <- sum(unlist(peaks$count))==0
}
sum(removePeaks)

peaks <- peaks[!removePeaks]
save(peaks,file=paste0(out, "_peaks.RData"))

########################
#### create indices ####
########################
peaks.df.index=data.frame(chr=seqnames(peaks),
                          start=start(peaks),
                          end=end(peaks),
                          width=width(peaks))
if(!is.null(score(peaks)))
  peaks.df.index$score=score(peaks)
# maximum height
peaks.df.index$height=unlist(lapply(peaks$count,max))
# area
peaks.df.index$area=unlist(lapply(peaks$count,sum))
# full width at half maximum
peaks.df.index$width.height.2=mapply(function(count,height){max(which(count>=height/2))-min(which(count>=height/2))},peaks$count,peaks.df.index$height)
# number of local peaks
if(!exists('N'))
  N=min(20,peaks.df.index$width) # distance knots
if(N>min(peaks.df.index$width)){
  warning('\'N\' is too big. Setting it to the minimum peak width.')
  N=min(peaks.df.index$width)
}
if(N<4){
  if(min(peaks.df.index$width)<4)
    stop('There are very short peaks (<4 nucleotides).')
  warning('\'N\' is too low. Setting to the default value.')
  N=min(20,peaks.df.index$width)
}
if(1 %in% unlist(lapply(peaks$count,function(count) length(unique(count)))))
  warning('There are peaks with constant coverage.')
if(!exists('toll'))
  toll <- min(50,peaks.df.index$width) # minimum distance between maxima
if(toll>min(peaks.df.index$width)){
  warning('\'toll\' is too big. Setting it to the minimum peak width.')
  toll=min(peaks.df.index$width)
}
toll2=0.2 # minimum high (in percentage) of maxima
peaks.df.index$local.peaks=get_local_peaks(peaks$count,N,toll,toll2)

# shape index M divided by the maximum height
peaks.df.index$M.height=unlist(lapply(peaks$count,function(count) get_M(length(count)+1,c(0,count))))/peaks.df.index$height

### SAVE INDICES
save(peaks.df.index,file=paste0(out, "_peaks_index.RData"))


# scatter plot
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    Cor <- abs(cor(x, y)) # Remove abs function if desired
    txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
    if(missing(cex.cor)) {
        cex.cor <- 0.4 / strwidth(txt)
    }
    text(0.5, 0.5, txt,
         cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}

## Plot indices correlation
pdf(paste0(out, "_scatterplot.pdf"),10,10)
pairs(peaks.df.index[,c("height","area","width.height.2","local.peaks","M.height")],
      col="black",main='Shape indices',pch=3,
      labels=expression(h,A,w [h/2],p[local],frac(M,h)),cex.labels=1.8, upper.panel = panel.cor, lower.panel = panel.smooth)
invisible(dev.off())

####################
#### clustering ####
####################


# standardization
index.std=peaks.df.index[,c("height","area","width.height.2","local.peaks","M.height")]

# Looking for NA
apply(index.std, 2, function(index){sum(is.na(index))})

# Change NA for 0
index.std$M.height[is.na(index.std$M.height)] <- 0

index.std=apply(index.std,2,function(index){(index-mean(index))/sd(index)})



# k-means
K=lapply(1:10,function(k){kmeans(index.std,k,iter.max=100,nstart=10)})
within=unlist(lapply(K,function(K){K$tot.withinss}))

## SAVE k-means
save(K,within,file=paste0(out, "_kmeans.RData"))


# total within sum of squares plot
pdf(paste0(out, "_tot_within_ss.pdf"),10,10)
plot(1:length(within),within/within[1],xlim=c(1,10),ylim=c(0,1),xlab='Number of clusters k',ylab='Total within SS (%)',main='K-means on shape indices',las=1,type='b')
axis(side=1,at=1:10)
invisible(dev.off())


####################
### Clustering regions
###################
peaksBed <- read.table(bed, sep = "\t")
vector <- K[[10]]$cluster
for(i in 1:10){
  peakClass <- peaksBed[vector==i,]
  path = paste0(out, "_peak_class",i,".bed")
  write.table(peakClass, path, sep = "\t", quote = F, row.names = F, col.names = F)
}


