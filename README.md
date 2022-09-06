# DESEq2_Analysis
####  Do different sample groups express genes differentially
#### reference 1:
#### https://github.com/nattoh03/RNASeq-Analysis-with-Hisat2/blob/main/pipeline_fastqc_alignHISAT2_bamsorted_output.txt

Three questions to answer

#### which genes do the reads (samples) belong to ?
    done in supercomputer, alignment of ref genome with Hisat2
    reference 1 above

#### How many reads align to a specifi gene?
   done in supercomputer,
   using htseq-count, 
   the gff file and 
   the sorted.bam file 
   to generate a txt/csv file for downstream analysis in R 
   reference 1 above
   
#### Do different samples express genes differentially?
    now this part is done in R
    using library(DESeq2)
    the .txt file from htseq-count
    and sample description file (.txt)
    
    
####To install this DESeq2 package, start R (version "4.2") and enter:

R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
####(source: https://bioconductor.org/packages/release/bioc/html/DESeq2.html)




#Load packages----
library(limma)
library(edgeR)
library(ggplot2)
library(tidyverse)
library(ggpubr)
#set the working directory. Go to compartment C on your computer and create a folder named Diana. And transfer the excel file into this folder, name the excel file as featurecount, And make sure this excel file is saved as comma delimited. 
getwd()
setwd("D:/new isaiah")
alldata <- read.csv("table_gene_counts.csv", header = T, row.names =1)
names(alldata)

allres_vs_susce <- alldata[, c(1, 3, 4, 5, 12, 13, 6, 7, 9, 11, 14)]
kombewa <- alldata[, c(1, 7, 6, 5, 4, 3, 2, 10)]
names(kombewa)

siaya <- alldata[, c(1, 9, 8, 10, 8, 2, 9)]
names(siaya)


port_vict <- alldata[, c(1, 11, 12, 2, 12, 10, 11)]
names(port_vict)


teso <- alldata[, c(1, 14, 13, 2, 13, 10, 14)]
names(teso)

metaTeso <- data.frame(samples =c(Teso_resistant),




library(DESeq2 )
library(ggplot2)

    
