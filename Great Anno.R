setwd("C:/Users/Alber/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu20.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/alberto/TFM/TGFb2h/SIC/results_ATAC")

library("rGREAT")



### Get the rGREAT jobs. 

  
job_ATAC_1 = submitGreatJob("TGFb2h_peak_class1.bed", species = "mm9", rule = "basalPlusExt", adv_upstream = 50, adv_downstream = 50, adv_span = 1000)
job_ATAC_2 = submitGreatJob("TGFb2h_peak_class2.bed", species = "mm9", rule = "basalPlusExt", adv_upstream = 50, adv_downstream = 50, adv_span = 1000)
job_ATAC_3 = submitGreatJob("TGFb2h_peak_class3.bed", species = "mm9", rule = "basalPlusExt", adv_upstream = 50, adv_downstream = 50, adv_span = 1000)
job_ATAC_4 = submitGreatJob("TGFb2h_peak_class4.bed", species = "mm9", rule = "basalPlusExt", adv_upstream = 50, adv_downstream = 50, adv_span = 1000)
job_ATAC_5 = submitGreatJob("TGFb2h_peak_class5.bed", species = "mm9", rule = "basalPlusExt", adv_upstream = 50, adv_downstream = 50, adv_span = 1000)
job_ATAC_6 = submitGreatJob("TGFb2h_peak_class6.bed", species = "mm9", rule = "basalPlusExt", adv_upstream = 50, adv_downstream = 50, adv_span = 1000)
job_ATAC_7 = submitGreatJob("TGFb2h_peak_class7.bed", species = "mm9", rule = "basalPlusExt", adv_upstream = 50, adv_downstream = 50, adv_span = 1000)
job_ATAC_8 = submitGreatJob("TGFb2h_peak_class8.bed", species = "mm9", rule = "basalPlusExt", adv_upstream = 50, adv_downstream = 50, adv_span = 1000)
job_ATAC_9 = submitGreatJob("TGFb2h_peak_class9.bed", species = "mm9", rule = "basalPlusExt", adv_upstream = 50, adv_downstream = 50, adv_span = 1000)
job_ATAC_10 = submitGreatJob("TGFb2h_peak_class10.bed", species = "mm9", rule = "basalPlusExt", adv_upstream = 50, adv_downstream = 50, adv_span = 1000)


# Check for the available ontologies for our data

availableOntologies(job_ATAC_1)
availableCategories(job_ATAC_1)

Tab_1 =getEnrichmentTables(job_ATAC_1, ontology = c("GO Biological Process"), download_by = "tsv")
Tab_2 =getEnrichmentTables(job_ATAC_2, ontology = c("GO Biological Process"), download_by = "tsv")
Tab_3 =getEnrichmentTables(job_ATAC_3, ontology = c("GO Biological Process"), download_by = "tsv")
Tab_4 =getEnrichmentTables(job_ATAC_4, ontology = c("GO Biological Process"), download_by = "tsv")
Tab_5 =getEnrichmentTables(job_ATAC_5, ontology = c("GO Biological Process"), download_by = "tsv")
Tab_6 =getEnrichmentTables(job_ATAC_6, ontology = c("GO Biological Process"), download_by = "tsv")
Tab_7 =getEnrichmentTables(job_ATAC_7, ontology = c("GO Biological Process"), download_by = "tsv")
Tab_8 =getEnrichmentTables(job_ATAC_8, ontology = c("GO Biological Process"), download_by = "tsv")
Tab_9 =getEnrichmentTables(job_ATAC_9, ontology = c("GO Biological Process"), download_by = "tsv")
Tab_10 =getEnrichmentTables(job_ATAC_10, ontology = c("GO Biological Process"), download_by = "tsv")


## Plot TSS distnce

plotRegionGeneAssociationGraphs(job_ATAC_1,2)
plotRegionGeneAssociationGraphs(job_ATAC_2,2)
plotRegionGeneAssociationGraphs(job_ATAC_3,2)
plotRegionGeneAssociationGraphs(job_ATAC_4,2)
plotRegionGeneAssociationGraphs(job_ATAC_5,2)
plotRegionGeneAssociationGraphs(job_ATAC_6,2)
plotRegionGeneAssociationGraphs(job_ATAC_7,2)
plotRegionGeneAssociationGraphs(job_ATAC_8,2)
plotRegionGeneAssociationGraphs(job_ATAC_9,2)
plotRegionGeneAssociationGraphs(job_ATAC_10,2)


setwd("C:/Users/Alber/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu20.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/alberto/TFM/TGFb2h/SIC/results_ATAC/Great Anno")

## Get the downloaded enriched biological processes from http://great.stanford.edu/public/html

ATAC_1_anno = read.table("BP_ATAC1.tsv", header = F, sep = "\t")
ATAC_2_anno = read.table("BP_ATAC2.tsv", header = F, sep = "\t")
ATAC_3_anno = read.table("BP_ATAC3.tsv", header = T, sep = "\t")
ATAC_4_anno = read.table("BP_ATAC4.tsv", header = F, sep = "\t")
ATAC_5_anno = read.table("BP_ATAC5.tsv", header = F, sep = "\t")
ATAC_6_anno = read.table("BP_ATAC6.tsv", header = F, sep = "\t")
ATAC_7_anno = read.table("BP_ATAC7.tsv", header = F, sep = "\t")
ATAC_8_anno = read.table("BP_ATAC8.tsv", header = F, sep = "\t")
ATAC_9_anno = read.table("BP_ATAC9.tsv", header = F, sep = "\t")
ATAC_10_anno = read.table("BP_ATAC10.tsv", header = F, sep = "\t")


colnames = c("Desc","BinomRank","Binom P-value","Binom FDR Q-value","Binom Fold Enrichment","Binom ObsRegions","Binom SetCov","HyperRank","Hyper FDR Q-value","Hyper Fold Enrichment","Hyper ObsRegions","Hyper total genes", "Hyper SetCov")
colnames(ATAC_10_anno)= colnames
write.table(ATAC_10_anno, file = "ATAC_peaks10_Biological_process.csv", sep = "\t")
save(ATAC_10_anno, file = "BP_anno_ATAC10.RData")
     

