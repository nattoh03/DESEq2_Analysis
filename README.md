# DESEq2_Analysis
####  Do different sample groups express genes differentially
#### reference 1:
#### https://github.com/nattoh03/RNASeq-Analysis-with-Hisat2/blob/main/pipeline_fastqc_alignHISAT2_bamsorted_output.txt

The Three questions to answer

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
    the .txt or .csv file from htseq-count
    and sample description file (.txt or .csv)
    
 
 
 
 
DOWNSTREAM ANALYSIS IN R 
 
#### To install this DESeq2 package, start R (version "4.2") and enter:
    R
    if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

    BiocManager::install("DESeq2")
#### (source: https://bioconductor.org/packages/release/bioc/html/DESeq2.html)


#### Load packages ----
    library(ggplot2)
    library(tidyverse)
    library(ggpubr)
    library(limma)
    library(DESeq2)

      
#### set the working directory. Go to compartment D on your computer and create a folder named isaiah. And transfer the csv file into this folder, name the excel file as table_gene_counts.csv, And make sure this excel file is saved as comma delimited. 
    # check your present working directory
    getwd()
    
    
    set a new folder as your working directory 
    
           setwd("D:/isaiah")
    
    
#### To load your data, and give it a new name e.g alldata
    
         alldata <- read.csv("table_gene_counts.csv", header = T, row.names =1)
           
#### To check the names of the columns in your data to guide you if you want to extract a part of it
    
         names(alldata)
    
    
    extract sample from one site e.g kombewa area (i.e. resistant vs susceptible)
   
    kombewa <- alldata[, c(1, 7, 6, 5, 4, 3, 2, 10)]
    names(kombewa)
    
    
#### To extract samples from siaya area  (i.e. resistant vs susceptible)
    
    siaya <- alldata[, c(1, 9, 8, 10, 8, 2, 9)]
    names(siaya)
    
    
    extract samples from port victoria area  (i.e. resistant vs susceptible)
    
    port_vict <- alldata[, c(1, 11, 12, 2, 12, 10, 11)]
    names(port_vict)


#### downstream analysis of Teso samples (resistant vs susceptible) # plotting of a volcano plot
    extract samples from  teso area  (i.e. resistant vs susceptible)
    
    teso <- alldata[, c(1, 14, 13, 2, 13, 10, 14)]
    names(teso)
    
    write.csv(teso, "teso_res_vs_susce.csv")
    
    teso_res_vs_susce <- as.matrix(read.csv("teso_res_vs_susce.csv", header = T, row.names = "gene_id"))
    
    teso_res_vs_susce <- teso_res_vs_susce[, c(2, 3, 4, 5, 6, 7)]
    
    sample info
    teso_res_vs_susce_Info2R2 <- as.matrix(read.delim("teso_res_vs_susceINFO.txt", header = T, sep = '\t', row.names =1))
   
   #### undertake differential analysis, but first ensure the library is loaded
    library(DESeq2)

    dds_teso<- DESeqDataSetFromMatrix(teso_res_vs_susce, teso_res_vs_susce_Info2R2, ~condition)
    dds_teso
    
    keep_teso <- rowSums(counts(dds_teso)) >20
    
    
    dds_teso <- dds_teso[keep_teso,]
    
    ddsDE_teso <- DESeq(dds_teso)
    
    
    results(ddsDE_teso, alpha = 0.05)
    
    normCounts <- counts(ddsDE_teso, normalized =T)
    
    write.csv(normCounts, "norma_teso_res_vs_sus.csv")
    
    
    Res_teso<- results(ddsDE_teso, alpha = 0.05)
    
    summary(Res_teso)
    
    resOrdered_all <- Res_teso[order(Res_all$padj),]
    
    write.csv(resOrdered_all, "deseq.teso_res_susce.csv")
    
    resultsNames(ddsDE_teso)
    
    plotMA(ddsDE_teso, ylim=c(-3,3))
    
    deseqRes_teso <- read.csv("deseq.teso_res_susce.csv", row.names = 1)
    
    deseqRes_teso <- na.omit(deseqRes_teso)
    
    deseqRes_teso$sig <- ifelse(deseqRes_teso$padj <=0.05, "yes", "no")
    
    library(ggplot2)
    
    ggplot(deseqRes_teso, aes(x=baseMean, y=log2FoldChange, col =sig)) + geom_point()
    
    ggplot(deseqRes_teso, aes(x=log2FoldChange, y=-log10(padj), col =sig)) + geom_point()
    
    ggplot(deseqRes_teso, aes(x=log2FoldChange, y=-log10(padj), col =sig)) + geom_point(stroke =2 , shape = 16, size=1)+ 
    theme(panel.background = element_rect(fill = "white", color = "grey50"))
    
    ggplot(deseqRes_teso, aes(x=log2FoldChange, y=-log10(padj), col =sig)) +
    geom_point(stroke =2 , shape = 16, size=1)+
    theme(panel.background = element_rect(fill = "white", color = "grey50")) +
    #scale_y_continuous(limits = c(0, 120))+
    geom_vline(xintercept=c(-1, 1), col="black", lty=3) + geom_hline(yintercept = -log10(0.01), lty=3, color="black")+
    theme(legend.title=element_blank())+
    #scale_y_continuous(limits = c(0, 90))+
    ggtitle("Teso-Res vs Susc")+
    labs(x = expression ("log"[2](FC)), y=expression("-log"[10](FDR)))
   
