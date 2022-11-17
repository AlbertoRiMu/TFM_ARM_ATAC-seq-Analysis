# TFM_ARM_ATAC-seq-
This project carried out the clustering and characterization of the data of ATAC-seq analysis. Using the ChIP-SIC software to cluster the reads based on their 
morphological characteristics.

After getting your downloaded sequences from ATAC-seq analysis you should process them with the https://github.com/kundajelab/atac_dnase_pipelines pipeline to get the .bam processed files
The general workflow of this project starts with the Create bw files.rmd to generate the 3 bw files of the inputs, ATAC, open chromatin regions and mononucleosome 
free regions.
Then the SIC-ChIP.R script will perform the clustering of the bw files, make sure to switch between the three types of input. 
The rest of the scripts will perform the characterization of the clusters, and annotate the biological information of each cluster, execute them in the following order.

- ViolinPlot.R
- Heatmap.R
- *PlotMatrix.R
- *Motifs.R
- GreatAnno..R

The last script being the logisticReg.R will create a prediction model based on the five indices to determine if the region is enhancer or promoter

*Please note that some of the commented code on some of the scripts will need to be run using Linux.
