## Get fimo files with motif list

#### GET FASTA FILES FROM .BED
## for i in TGFb2h_peak_class10  TGFb2h_peak_class1  TGFb2h_peak_class2  TGFb2h_peak_class3  TGFb2h_peak_class4  TGFb2h_peak_class5  TGFb2h_peak_class6  TGFb2h_peak_class7  TGFb2h_peak_class8  TGFb2h_peak_class9; do 
## bedtools getfasta -fi ../motivos/SIC/mm9_filtered.fa -bed $i.bed -fo $i.fa & 
## done

## GET BACKGROUND SEQUENCES
## for i in TGFb2h_peak_class10  TGFb2h_peak_class1  TGFb2h_peak_class2  TGFb2h_peak_class3  TGFb2h_peak_class4  TGFb2h_peak_class5  TGFb2h_peak_class6  TGFb2h_peak_class7  TGFb2h_peak_class8  TGFb2h_peak_class9; do scrambleFasta.pl $i.fa > $i.bg.fa & done

## for i in TGFb2h_peak_class10  TGFb2h_peak_class1  TGFb2h_peak_class2  TGFb2h_peak_class3  TGFb2h_peak_class4  TGFb2h_peak_class5  TGFb2h_peak_class6  TGFb2h_peak_class7  TGFb2h_peak_class8  TGFb2h_peak_class9; do 
##fimo --norc --max-stored-scores 100000000 -o ./${i}_mot ../motivos/SIC/MOTIF_MEME.txt ${i}.fa &
##done

## for i in TGFb2h_peak_class10  TGFb2h_peak_class1  TGFb2h_peak_class2  TGFb2h_peak_class3  TGFb2h_peak_class4  TGFb2h_peak_class5  TGFb2h_peak_class6  TGFb2h_peak_class7  TGFb2h_peak_class8  TGFb2h_peak_class9; do 
##	fimo --norc --max-stored-scores 100000000 -o ./${i}_bg ../motivos/SIC/MOTIF_MEME.txt ${i}.bg.fa & 
##	done

## rm */*.gff */*.xml */*.html





setwd("Escritorio/SIC/results")



### Create dirs for fasta and bg sequences

dirs <-  dir(pattern = "_mot")
dirs_bg <- dir(pattern = "_bg")

names <- gsub(pattern = "TGFb2h_peak_", replacement = "", x = dirs)
names_bg <- gsub(pattern = "TGFb2h_peak_", replacement = "", x = dirs_bg)

### Read tables

for(i in 1:length(names)){
  assign(names[i], read.table(paste0(dirs[i], "/fimo.tsv"), header = T))
  assign(names_bg[i], read.table(paste0(dirs_bg[i], "/fimo.tsv"), header = T))
}


#### Filter sequences

cutoff <- 12
filtScore <- function(tab){
  return(tab[tab$score > cutoff,])
}

names_filt <- paste0(names, "_filt")
names_bg_filt <- paste0(names_bg, "_filt")


for(i in 1:length(names)){
  assign(names_filt[i], filtScore(get(names[i])))
  assign(names_bg_filt[i], filtScore(get(names_bg[i])))
}


#### Frequency tables


names_freq <- paste0(names, "_freq")
names_bg_freq <- paste0(names_bg, "_freq")

for(i in 1:length(names)){
  assign(names_freq[i], as.data.frame(table(get(names_filt[i])$motif_alt_id)))
  assign(names_bg_freq[i], as.data.frame(table(get(names_bg_filt[i])$motif_alt_id)))
}


######### ENRICHMENT

motifs <- unique(as.character(unlist(lapply(names_freq, function(x) {get(x)[,1]})),
                 unlist(lapply(names_bg_freq, function(x) {get(x)[,1]}))))

lookMot <- function(freq, motif){
  freq[,1] <- as.character(freq[,1])
  idx <- match(motif, freq[,1])
  if(!is.na(idx)) {
    return(as.numeric(freq[idx,2]))
  } else {
    return(1)
  }
}


names_motifs <- paste0(names, "_motifs")
names_bg_motifs <- paste0(names_bg, "_motifs")

# Get motifs for fasta seq
v <- c()
for(i in 1:length(names)){
  tab <- get(names_freq[i])
  for(j in 1:length(motifs)){
    v[j] <- lookMot(tab, motifs[j])
  }
  names(v) <- motifs
  assign(names_motifs[i], v)
}
# Get motifs for bg seq
v <- c()
for(i in 1:length(names)){
  tab <- get(names_bg_freq[i])
  for(j in 1:length(motifs)){
    v[j] <- lookMot(tab, motifs[j])
  }
  names(v) <- motifs
  assign(names_bg_motifs[i], v)
}


names_enr <- paste0(names, "_enr")

for(i in 1:length(names)){
  class  <- get(names_motifs[i])+5
  class_bg <- get(names_bg_motifs[i])+5
  neg <-  data.frame()
  neg_bg <- data.frame()
  for(j in c(1:length(names))[-i]){
    neg <- rbind(neg, get(names_motifs[j]))
    neg_bg <- rbind(neg_bg, get(names_bg_motifs[j]))
  }
  neg <- colSums(neg)+5
  neg_bg <- colSums(neg_bg)+5
  enr <- (class/(class_bg/5))/(neg/(neg_bg/5)) # divide by 5 because the bg sequences were 5 times more
  assign(names_enr[i], enr)
}


### P-VALUE

names_p <- paste0(names, "_p")

for(i in 1:length(names)){
  class  <- get(names_motifs[i])
  class_bg <- get(names_bg_motifs[i])
  neg <-  data.frame()
  neg_bg <- data.frame()
  for(j in c(1:length(names))[-i]){
    neg <- rbind(neg, get(names_motifs[j]))
    neg_bg <- rbind(neg_bg, get(names_bg_motifs[j]))
  }
  neg <- colSums(neg)
  neg_bg <- colSums(neg_bg)
  p <- c()
  for(m in 1:length(class)){
    p[m] <- fisher.test(matrix(c(class[m], class_bg[m], neg[m], neg_bg[m]), 
                               nrow = 2, ncol = 2))$p.value
  }
  assign(names_p[i], p)
}

# Create enriched motifs table with p-value and enrichment value
enr_tab <- data.frame()
p_tab <- data.frame()
for(i in 1:length(names)){
  enr_tab <- rbind(enr_tab, get(names_enr[i]))
  p_tab <- rbind(p_tab, get(names_p[i]))
}

enr_tab <- t(enr_tab)
rownames(enr_tab) <- motifs
p_tab <- t(p_tab)
rownames(p_tab) <- motifs


#### FILTER BY P-VALUE AND ENRICHMENT

names_mot_filt <- paste0(names, "_mot_filt")

filtrarMotivos <- function(enrich, pval){
  tabla <- cbind(enrich, pval)
  tabla2 <- tabla[enrich>1.5 & pval<0.05,]
  return(tabla2)
}

for(i in 1:length(names)){
  assign(names_mot_filt[i], filtrarMotivos(get(names_enr[i]), get(names_p[i])))
}

get(names_mot_filt[6])


## PLOT
library(ggplot2)
library(ggrepel)
library(ggpubr)

plotMotifs <- function(tab, title){
  tab <- as.data.frame(tab)
  tab <- tab[order(tab$enrich, decreasing = T),]
  if(nrow(tab) > 20){
    tab <- tab[c(1:20),]
  }
  .labs <- rownames(tab)
  Mplot <- ggplot(tab, aes(x = tab[,1], y = -log10(tab[,2])))
  Mplot <- Mplot + geom_point() +
    geom_text_repel(aes(label = .labs), size = 5) + theme_bw() +
    ggtitle(title) +
    labs(y = "-log10(p-value)", x = "Enrichment") + theme(plot.title = element_text(hjust = 0.5)) 
  return(Mplot)
}

p01 = plotMotifs(get(names_mot_filt[1]),gsub(pattern = "_mot",replacement = "",names[1]))
p02 = plotMotifs(get(names_mot_filt[2]),gsub(pattern = "_mot",replacement = "",names[2]))
p03 = plotMotifs(get(names_mot_filt[3]),gsub(pattern = "_mot",replacement = "",names[3]))
p04 = plotMotifs(get(names_mot_filt[4]),gsub(pattern = "_mot",replacement = "",names[4]))
p05 = plotMotifs(get(names_mot_filt[5]),gsub(pattern = "_mot",replacement = "",names[5]))
p06 = plotMotifs(get(names_mot_filt[6]),gsub(pattern = "_mot",replacement = "",names[6]))
p07 = plotMotifs(get(names_mot_filt[7]),gsub(pattern = "_mot",replacement = "",names[7]))
p08 = plotMotifs(get(names_mot_filt[8]),gsub(pattern = "_mot",replacement = "",names[8]))
p09 = plotMotifs(get(names_mot_filt[9]),gsub(pattern = "_mot",replacement = "",names[9]))
p10 = plotMotifs(get(names_mot_filt[10]),gsub(pattern = "_mot",replacement = "",names[10]))




plot <- ggarrange(p02, p03, p04, p05,p06, p07, p08, p10, nrow = 2, ncol = 5)


pdf("Motif_enrichment_promoters.pdf")
plot
dev.off()


