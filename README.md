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
    
    
#### To install this DESeq2 package, start R (version "4.2") and enter:
    R
    if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

    BiocManager::install("DESeq2")
#### (source: https://bioconductor.org/packages/release/bioc/html/DESeq2.html)




#### Load packages----
#### 
    library(edgeR)
    library(ggplot2)
    library(tidyverse)
    library(ggpubr)
    library(limma)
    
    


#### set the working directory. Go to compartment D on your computer and create a folder named isaiah. And transfer the csv file into this folder, name the excel file as table_gene_counts.csv, And make sure this excel file is saved as comma delimited. 
    getwd()
    library(edgeR)
    library(ggplot2)
    setwd("D:/isaiah")
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

    
############## ALL RES   VS   ALL SUSCE 
allres_vs_susce1 <- alldata[, c(1, 3, 4, 5, 12, 13, 6, 7, 9, 11, 14)]

names(allres_vs_susce1)

write.csv(allres_vs_susce1, "allres_vs_susce.csv")

allres_vs_susce <- as.matrix(read.csv("allres_vs_susce.csv", header = T, row.names = "gene_id"))
allres_vs_susce <- allres_vs_susce[, c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11)]
names(allres_vs_susce)

allres_vs_susce_Info2R2 <- as.matrix(read.delim("allresvssus.txt", header = T, sep = '\t', row.names =1))
names(allres_vs_susce_Info2R2)

library(DESeq2)

dds_all<- DESeqDataSetFromMatrix(allres_vs_susce, allres_vs_susce_Info2R2, ~condition)
dds_all
### if gene wide estimation dispersio are within 2 order
## dds <- estimateDispersionsGeneEst(dds)
# dispersions(dds) <- mcols(dds)$dispGeneEst
#...then continue with testing using nbinomWaldTest or nbinomLRT

keep_all <- rowSums(counts(dds_all)) >50
dds_all <- dds_all[keep_all,]

ddsDE_all <- DESeq(dds_all)
results(ddsDE_all, alpha = 0.05)

normCounts <- counts(ddsDE_all, normalized =T)
write.csv(normCounts, "norma_all_RES_SUSC.csv")

Res_all<- results(ddsDE_all, alpha = 0.05)
summary(Res_all)
resOrdered_all <- Res_all[order(Res_all$padj),]
write.csv(resOrdered_all, "deseq.all_res_susce.csv")
resultsNames(ddsDE_all)
plotMA(ddsDE_all, ylim=c(-3,3))


deseqRes_all <- read.csv("deseq.all_res_susce.csv", row.names = 1)
deseqRes_all <- na.omit(deseqRes_all)
deseqRes_all$sig <- ifelse(deseqRes_all$padj <=0.05, "yes", "no")
library(ggplot2)
ggplot(deseqRes_all, aes(x=baseMean, y=log2FoldChange, col =sig)) +
  geom_point()

ggplot(deseqRes_all, aes(x=log2FoldChange, y=-log10(padj), col =sig)) +
  geom_point()


ggplot(deseqRes_all, aes(x=log2FoldChange, y=-log10(padj), col =sig)) +
  geom_point(stroke =2 , shape = 16, size=1)+
  theme(panel.background = element_rect(fill = "white", color = "grey50"))



ggplot(deseqRes_all, aes(x=log2FoldChange, y=-log10(padj), col =sig)) +
  geom_point(stroke =2 , shape = 16, size=1)+
  theme(panel.background = element_rect(fill = "white", color = "grey50")) +
  #scale_y_continuous(limits = c(0, 120))+
  geom_vline(xintercept=c(-1, 1), col="black", lty=3) + geom_hline(yintercept = -log10(0.01), lty=3, color="black")+
  theme(legend.title=element_blank())+
  #scale_y_continuous(limits = c(0, 90))+
  ggtitle("I-MAL vs Rock")+
  labs(x = expression ("log"[2](FC)), y=expression("-log"[10](FDR)))
