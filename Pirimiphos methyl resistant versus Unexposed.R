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

#### 	Pirimiphos -methyl resistant versus Unexposed/control (PM 
pm_res_vs_unexp <- a_data[, c(1, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75)]
names(pm_res_vs_unexp)
write.csv(pm_res_vs_unexp, "Pirimiphos_resistant_vs_unexposed.csv")


head(pm_res_vs_unexp)
rownames(pm_res_vs_unexp) <- pm_res_vs_unexp[,1]
pm_res_vs_unexp <- pm_res_vs_unexp[,-1]
head(sus_kis_pm_res)

all(colnames(pm_res_vs_unexp) == matrix_Info2R2$sample.id)

str(pm_res_vs_unexp)


#sample info
matrix_Info2R2 <- as.matrix(read.delim("pirimiphos_unexposed_info.txt", header = T, sep = '\t', row.names =1))

all(colnames(pm_res_vs_unexp) %in% rownames(matrix_Info2R2))

all(colnames(pm_res_vs_unexp) == rownames(matrix_Info2R2))


#undertake differential analysis, but first ensure the library is loaded
library(DESeq2)

dd_pm_res_vs_unexp<- DESeqDataSetFromMatrix(pm_res_vs_unexp, matrix_Info2R2, ~condition)
dd_pm_res_vs_unexp

keep_pm_res_vs_unexp <- rowSums(counts(dd_pm_res_vs_unexp)) >20


dd_pm_res_vs_unexp <- dd_pm_res_vs_unexp[keep_pm_res_vs_unexp]

ddsDE_pm_res_vs_unexp <- DESeq(dd_pm_res_vs_unexp)


results(ddsDE_pm_res_vs_unexp, alpha = 0.05)

normCounts <- counts(ddsDE_pm_res_vs_unexp, normalized =T)

write.csv(normCounts, "normCount_pirimiphos_res_vs_unexposed.csv")


Res_pm_res_vs_unexp<- results(ddsDE_pm_res_vs_unexp, alpha = 0.05)

summary(Res_pm_res_vs_unexp)

resOrdered_all <- Res_pm_res_vs_unexp[order(Res_pm_res_vs_unexp$padj),]

write.csv(resOrdered_all, "deseq.pirimiphos_res_vs_unexposed.csv")

resultsNames(ddsDE_pm_res_vs_unexp)

plotMA(ddsDE_pm_res_vs_unexp, ylim=c(-3,3))

deseqRes_pm_res_vs_unexp <- read.csv("deseq.pirimiphos_res_vs_unexposed.csv", row.names = 1)

deseqRes_pm_res_vs_unexp <- na.omit(deseqRes_pm_res_vs_unexp)

deseqRes_pm_res_vs_unexp$differentialexpressed <- 'No'

deseqRes_pm_res_vs_unexp$differentialexpressed[deseqRes_pm_res_vs_unexp$log2FoldChange>0.06&deseqRes_pm_res_vs_unexp$padj<0.05] <- 'Up'
deseqRes_pm_res_vs_unexp$differentialexpressed[deseqRes_pm_res_vs_unexp$log2FoldChange<0.06&deseqRes_pm_res_vs_unexp$padj<0.05] <- 'Down'
####
deseqRes_sus_kis_perm_re$differentialexpressed[deseqRes_sus_kis_perm_re$log2FoldChange>0.6&deseqRes_sus_kis_perm_re$padj<0.05] <- 'Up'
deseqRes_sus_kis_perm_re$differentialexpressed[deseqRes_sus_kis_perm_re$log2FoldChange<0.6&deseqRes_sus_kis_perm_re$padj<0.05] <- 'Down'



### deseqRes_sus_kis_perm_re$sig <- ifelse(deseqRes_sus_kis_perm_re$padj <=0.05, "yes", "no")

library(ggplot2)


ggplot(deseqRes_pm_res_vs_unexp, aes(x=log2FoldChange, y=-log10(padj), col =differentialexpressed)) + geom_point()

ggplot(deseqRes_pm_res_vs_unexp, aes(x=log2FoldChange, y=-log10(padj), col =differentialexpressed)) + geom_point(stroke =2 , shape = 16, size=1)+ 
  theme(panel.background = element_rect(fill = "white", color = "grey50"))

ggplot(deseqRes_pm_res_vs_unexp, aes(x=log2FoldChange, y=-log10(padj), col =differentialexpressed)) +
  geom_point(stroke =2 , shape = 16, size=1)+
  theme(panel.background = element_rect(fill = "white", color = "grey50")) +
  #scale_y_continuous(limits = c(0, 120))+
  geom_vline(xintercept=c(-1, 1), col="black", lty=3) + geom_hline(yintercept = -log10(0.01), lty=3, color="black")+
  theme(legend.title=element_blank())+
  #scale_y_continuous(limits = c(0, 90))+
  ggtitle("Pirimiphos resistant vs unexposed")+
  labs(x = expression ("log"[2](FC)), y=expression("-log"[10](FDR)))




#####

p1 <- ggplot(object, aes(CTRL, CIS)) + geom_point() + ggtitle("Macrophage CTRL vs CIS") + theme(plot.title = element_text(hjust = 0.5)) + geom_point(size = 1,alpha = 0.6) +   
  geom_smooth(method=lm, se=TRUE, color="brown", linetype="solid", size=1) 

p1 <- LabelPoints(plot = p1, points = c(genes.to.labelCtrl, genes.to.labelCis), color="blue", repel = TRUE, xnudge=0, ynudge=0)

plot_grid(p1)




