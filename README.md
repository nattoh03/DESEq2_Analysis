# DESEq2_Analysis

####  Do different sample groups express genes differentially
#### reference 1:
## The Three questions to answer

#### which genes do the reads (samples) belongs to ?
    done in supercomputer, alignment of ref genome with Hisat2
    

#### How many reads align to a specific gene?
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
 ### in R-studio
 
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

      
#### set the working directory. 
Go to compartment D on your computer and create a folder named isaiah. And transfer the csv file into this folder, name the excel file as table_gene_counts.csv, And make sure this excel file is saved as comma delimited.
#### check your present working directory
getwd()

#### or set a new folder as your working directory using this command

setwd("C:/Users/nattohz/E/RNASeq")


    
    
#### To load your data, and give it a new name e.g alldata
    
         alldata <- read.csv("table_gene_counts.csv", header = T, row.names =1)
           
#### To check the names of the columns in your data to guide you if you want to extract a part of it
    
         names(alldata)
    
    
   
 #To load your data, and give it a new name e.g alldata
a_data <- read.csv("data_feature_count.csv", header = T, row.names =1)
names(a_data)
#### To extract specific samples from a specific site or trated with a given insecticide  e.g perm_delta_pirim_vs_kis area (i.e. all resistant vs susceptible)

perm_delta_pirim_vs_kis <- a_data[, c(1, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 2, 3, 4, 5, 6, 7, 8, 9, 10, 21, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)]
names(perm_delta_pirim_vs_kis)

write.csv(perm_delta_pirim_vs_kis, "permethrin_vs_deltamethrin_vs_pirimi_vs_susceptible_kisumu.csv")

perm_delta_pirim_vs_kis

### check names 
head(perm_delta_pirim_vs_kis)
rownames(perm_delta_pirim_vs_kis) <- perm_delta_pirim_vs_kis[,1]
perm_delta_pirim_vs_kis <- perm_delta_pirim_vs_kis[,-1]
head(perm_delta_pirim_vs_kis)

#### check column names
all(colnames(perm_delta_pirim_vs_kis) == matrix_Info2R2$sample.id)

str(perm_delta_pirim_vs_kis)


#### sample info
matrix_Info2R2 <- as.matrix(read.delim("perm_delta_piri_vs_kisumu_info.txt", header = T, sep = '\t', row.names =1))
#

all(colnames(perm_delta_pirim_vs_kis) %in% rownames(matrix_Info2R2))
#

all(colnames(perm_delta_pirim_vs_kis) == rownames(matrix_Info2R2))


#### undertake differential analysis, but first ensure the library is loaded
library(DESeq2)

dd_perm_delta_pirim_vs_kis<- DESeqDataSetFromMatrix(perm_delta_pirim_vs_kis, matrix_Info2R2, ~condition)
dd_perm_delta_pirim_vs_kis

keep_dd_perm_delta_pirim_vs_kis <- rowSums(counts(dd_perm_delta_pirim_vs_kis)) >20


dd_perm_delta_pirim_vs_kis <- dd_sus_kis_perm_res[keep_dd_perm_delta_pirim_vs_kis]

ddsDE_perm_delta_pirim_vs_kis <- DESeq(dd_perm_delta_pirim_vs_kis)


results(ddsDE_perm_delta_pirim_vs_kis, alpha = 0.05)

normCounts <- counts(ddsDE_perm_delta_pirim_vs_kis, normalized =T)

write.csv(normCounts, "norma_perm_delta_pirim_vs_kis.csv")


Res_perm_delta_pirim_vs_kis<- results(ddsDE_perm_delta_pirim_vs_kis, alpha = 0.05)

summary(Res_perm_delta_pirim_vs_kis)

resOrdered_all <- Res_perm_delta_pirim_vs_kis[order(Res_perm_delta_pirim_vs_kis$padj),]

write.csv(resOrdered_all, "deseq.perm_delta_pirim_vs_kis.csv")

resultsNames(ddsDE_perm_delta_pirim_vs_kis)

plotMA(ddsDE_perm_delta_pirim_vs_kis, ylim=c(-3,3))

deseqRes_perm_delta_pirim_vs_kis <- read.csv("deseq.perm_delta_pirim_vs_kis.csv", row.names = 1)

deseqRes_perm_delta_pirim_vs_kis <- na.omit(deseqRes_perm_delta_pirim_vs_kis)

deseqRes_perm_delta_pirim_vs_kis$differentialexpressed <- 'No'

deseqRes_perm_delta_pirim_vs_kis$differentialexpressed[deseqRes_perm_delta_pirim_vs_kis$log2FoldChange>0.06&deseqRes_perm_delta_pirim_vs_kis$padj<0.05] <- 'Up'
deseqRes_perm_delta_pirim_vs_kis$differentialexpressed[deseqRes_perm_delta_pirim_vs_kis$log2FoldChange<0.06&deseqRes_perm_delta_pirim_vs_kis$padj<0.05] <- 'Down'
write.csv(deseqRes_perm_delta_pirim_vs_kis, "diff_deseq.perm_delta_pirim_vs_kis.csv")

####

ggplot(deseqRes_perm_delta_pirim_vs_kis, aes(x=log2FoldChange, y=-log10(padj), col =differentialexpressed)) + geom_point()

ggplot(deseqRes_perm_delta_pirim_vs_kis, aes(x=log2FoldChange, y=-log10(padj), col =differentialexpressed)) + geom_point(stroke =2 , shape = 16, size=1)+ 
  theme(panel.background = element_rect(fill = "white", color = "grey50"))

ggplot(deseqRes_perm_delta_pirim_vs_kis, aes(x=log2FoldChange, y=-log10(padj), col =differentialexpressed)) +
  geom_point(stroke =2 , shape = 16, size=1)+
  theme(panel.background = element_rect(fill = "white", color = "grey50")) +
  #scale_y_continuous(limits = c(0, 120))+
  geom_vline(xintercept=c(-1, 1), col="black", lty=3) + geom_hline(yintercept = -log10(0.01), lty=3, color="black")+
  theme(legend.title=element_blank())+
  #scale_y_continuous(limits = c(0, 90))+
  ggtitle("Perm vs Deltam vs Pirim vs Kisumu Susceptible")+
  labs(x = expression ("log"[2](FC)), y=expression("-log"[10](FDR)))
   
