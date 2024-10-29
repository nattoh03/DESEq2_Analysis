
#### kisumu vs unexposed


#Load packages ----
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(limma)
library(DESeq2)
install_github("cran/limma")



#set the working directory. Go to compartment D on your computer and create a folder named isaiah. And transfer the csv file into this folder, name the excel file as table_gene_counts.csv, And make sure this excel file is saved as comma delimited.
# check your present working directory
getwd()

#or set a new folder as your working directory 

setwd("C:/Users/indindag/OneDrive - State of Connecticut/Desktop/Godfrey/Isaiah 2nd Ms")
#To load your data, and give it a new name e.g alldata
a_data <- read.csv("data_feature_count.csv", header = T, row.names =1)
names(a_data)
unexposed_kis <- a_data[, c(1, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 2, 3, 4, 5, 6, 7, 8, 9, 10, 21, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50)]
names(unexposed_kis)
write.csv(unexposed_kis, "unexposed_vs_susceptible_kisumu.csv")



head(unexposed_kis)
rownames(unexposed_kis) <- unexposed_kis[,1]
unexposed_kis <- unexposed_kis[,-1]
head(unexposed_kis)

all(colnames(unexposed_kis) == matrix_Info2R2$sample.id)

str(unexposed_kis)


#sample info
matrix_Info2R2 <- as.matrix(read.delim("kis_susc_unexposed.txt", header = T, sep = '\t', row.names =1))

all(colnames(unexposed_kis) %in% rownames(matrix_Info2R2))

all(colnames(unexposed_kis) == rownames(matrix_Info2R2))


#undertake differential analysis, but first ensure the library is loaded
library(DESeq2)

dd_unexposed_kis<- DESeqDataSetFromMatrix(unexposed_kis, matrix_Info2R2, ~condition)
dd_unexposed_kis

keep_dd_unexposed_kis <- rowSums(counts(dd_unexposed_kis)) >20


dd_unexposed_kis <- dd_unexposed_kis[keep_dd_unexposed_kis]

ddsDE_unexposed_kis <- DESeq(dd_unexposed_kis)


results(ddsDE_unexposed_kis, alpha = 0.05)

normCounts <- counts(ddsDE_unexposed_kis, normalized =T)

write.csv(normCounts, "norma_kisumu_susce_vs_unexposed.csv")


Res_unexposed_kis<- results(ddsDE_unexposed_kis, alpha = 0.05)

summary(Res_unexposed_kis)

resOrdered_all <- Res_unexposed_kis[order(Res_unexposed_kis$padj),]

write.csv(resOrdered_all, "deseq.sus_kis_unexposed.csv")

resultsNames(ddsDE_unexposed_kis)

plotMA(ddsDE_unexposed_kis, ylim=c(-3,3))

deseqRes_unexposed_kis <- read.csv("deseq.sus_kis_unexposed.csv", row.names = 1)

deseqRes_unexposed_kis <- na.omit(deseqRes_unexposed_kis)

deseqRes_unexposed_kis$differentialexpressed <- 'No'

deseqRes_unexposed_kis$differentialexpressed[deseqRes_unexposed_kis$log2FoldChange>0.06&deseqRes_unexposed_kis$padj<0.05] <- 'Up'
deseqRes_unexposed_kis$differentialexpressed[deseqRes_unexposed_kis$log2FoldChange<0.06&deseqRes_unexposed_kis$padj<0.05] <- 'Down'
####



library(ggplot2)


ggplot(deseqRes_unexposed_kis, aes(x=log2FoldChange, y=-log10(padj), col =differentialexpressed)) + geom_point()

ggplot(deseqRes_unexposed_kis, aes(x=log2FoldChange, y=-log10(padj), col =differentialexpressed)) + geom_point(stroke =2 , shape = 16, size=1)+ 
  theme(panel.background = element_rect(fill = "white", color = "grey50"))

ggplot(deseqRes_unexposed_kis, aes(x=log2FoldChange, y=-log10(padj), col =differentialexpressed)) +
  geom_point(stroke =2 , shape = 16, size=1)+
  theme(panel.background = element_rect(fill = "white", color = "grey50")) +
  #scale_y_continuous(limits = c(0, 120))+
  geom_vline(xintercept=c(-1, 1), col="black", lty=3) + geom_hline(yintercept = -log10(0.01), lty=3, color="black")+
  theme(legend.title=element_blank())+
  #scale_y_continuous(limits = c(0, 90))+
  ggtitle("Kisumu strain vs Unexposed")+
  labs(x = expression ("log"[2](FC)), y=expression("-log"[10](FDR)))
