###Pirimiphos -methyl resistant versus Susceptible Kisumu Strain (PM vs KS)

#Load packages ----
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(limma)
library(DESeq2)
install_github("cran/limma")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

ABiocManager::install("limma")

#set the working directory. Go to compartment D on your computer and create a folder named isaiah. And transfer the csv file into this folder, name the excel file as table_gene_counts.csv, And make sure this excel file is saved as comma delimited.
# check your present working directory
getwd()

#or set a new folder as your working directory 

setwd("C:/Users/indindag/OneDrive - State of Connecticut/Desktop/Godfrey/Isaiah 2nd Ms")
#To load your data, and give it a new name e.g alldata
a_data <- read.csv("data_feature_count.csv", header = T, row.names =1)
names(a_data)
pirimo_res_vs_kis <- a_data[, c(1, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 2, 3, 4, 5, 6, 7, 8, 9, 10, 21, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)]
names(pirimo_res_vs_kis)
write.csv(pirimo_res_vs_kis, "pirimiphos_res_vs_kisumu.csv")

head(pirimo_res_vs_kis)
delta_kis_1 <- a_data[, c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48)]
names(delta_kis_1)


head(pirimo_res_vs_kis)
rownames(pirimo_res_vs_kis) <- pirimo_res_vs_kis[,1]
pirimo_res_vs_kis <- pirimo_res_vs_kis[,-1]
head(pirimo_res_vs_kis)

all(colnames(sus_kis_pm_res) == matrix_Info2R2$sample.id)

str(pirimo_res_vs_kis)


#sample info
matrix_Info2R2 <- as.matrix(read.delim("pirimiphos_kis_info.txt", header = T, sep = '\t', row.names =1))

all(colnames(pirimo_res_vs_kis) %in% rownames(matrix_Info2R2))

all(colnames(pirimo_res_vs_kis) == rownames(matrix_Info2R2))


#undertake differential analysis, but first ensure the library is loaded
library(DESeq2)

dd_pirimo_res_vs_kis<- DESeqDataSetFromMatrix(pirimo_res_vs_kis, matrix_Info2R2, ~condition)
dd_pirimo_res_vs_kis

keep_pirimo_res_vs_kis <- rowSums(counts(dd_pirimo_res_vs_kis)) >20


dd_pirimo_res_vs_kis <- dd_pirimo_res_vs_kis[keep_pirimo_res_vs_kis]

ddsDE_pirimo_res_vs_kis <- DESeq(dd_pirimo_res_vs_kis)


results(ddsDE_pirimo_res_vs_kis, alpha = 0.05)

normCounts <- counts(ddsDE_pirimo_res_vs_kis, normalized =T)

write.csv(normCounts, "norma_pirimiphos_vs_kisumu_normaCount.csv")


Res_pirimo_res_vs_kis<- results(ddsDE_pirimo_res_vs_kis, alpha = 0.05)

summary(Res_pirimo_res_vs_kis)

resOrdered_all_pirimo_res_vs_kis <- Res_pirimo_res_vs_kis[order(Res_pirimo_res_vs_kis$padj),]

write.csv(resOrdered_all_pirimo_res_vs_kis, "deseq.pirimiphos_kis.csv")

resultsNames(ddsDE_pirimo_res_vs_kis)

plotMA(ddsDE_pirimo_res_vs_kis, ylim=c(-3,3))

deseqRes_pirimo_res_vs_kis <- read.csv("deseq.pirimiphos_kis.csv", row.names = 1)

deseqRes_pirimo_res_vs_kis <- na.omit(deseqRes_pirimo_res_vs_kis)

deseqRes_pirimo_res_vs_kis$differentialexpressed <- 'No'

deseqRes_pirimo_res_vs_kis$differentialexpressed[deseqRes_pirimo_res_vs_kis$log2FoldChange>0.06&deseqRes_pirimo_res_vs_kis$padj<0.05] <- 'Up'
deseqRes_pirimo_res_vs_kis$differentialexpressed[deseqRes_pirimo_res_vs_kis$log2FoldChange<0.06&deseqRes_pirimo_res_vs_kis$padj<0.05] <- 'Down'
####
deseqRes_delta_kis$differentialexpressed[deseqRes_delta_kis$log2FoldChange>0.6&deseqRes_delta_kis$padj<0.05] <- 'Up'
deseqRes_delta_kis$differentialexpressed[deseqRes_delta_kis$log2FoldChange<0.6&deseqRes_delta_kis$padj<0.05] <- 'Down'



deseqRes_sus_kis_perm_re$sig <- ifelse(deseqRes_sus_kis_perm_re$padj <=0.05, "yes", "no")

library(ggplot2)

ggplot(deseqRes_pirimo_res_vs_kis, aes(x=baseMean, y=log2FoldChange, col =sig)) + geom_point()

ggplot(deseqRes_pirimo_res_vs_kis, aes(x=log2FoldChange, y=-log10(padj), col =differentialexpressed)) + geom_point()

ggplot(deseqRes_pirimo_res_vs_kis, aes(x=log2FoldChange, y=-log10(padj), col =differentialexpressed)) + geom_point(stroke =2 , shape = 16, size=1)+ 
  theme(panel.background = element_rect(fill = "white", color = "grey50"))

ggplot(deseqRes_pirimo_res_vs_kis, aes(x=log2FoldChange, y=-log10(padj), col =differentialexpressed)) +
  geom_point(stroke =2 , shape = 16, size=1)+
  theme(panel.background = element_rect(fill = "white", color = "grey50")) +
  #scale_y_continuous(limits = c(0, 120))+
  geom_vline(xintercept=c(-1, 1), col="black", lty=3) + geom_hline(yintercept = -log10(0.01), lty=3, color="black")+
  theme(legend.title=element_blank())+
  #scale_y_continuous(limits = c(0, 90))+
  ggtitle("Pirimiphos Resistant vs Kisumu Susceptible")+
  labs(x = expression ("log"[2](FC)), y=expression("-log"[10](FDR)))





