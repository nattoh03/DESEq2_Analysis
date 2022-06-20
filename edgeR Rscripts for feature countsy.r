#Load package-----
library(limma) # differential gene expression analysis
library(edgeR) # differential gene expression analysis
library(ggplot2) # for regular plotting
library(ggpubr) #  to plot figures in one window
library(data.table) # good to save big data frame in a file

##***-----SD vs DO-----
#load the data----
gp <- read.csv("annotation.csv")
x<-read.delim("Aarab.antisense.count_dd.txt",skip=0, sep="\t", check.names=FALSE, row.names = 1)

#indexing the columns of interest
counts <- x[,c('D01',	'D02','D03', "SD01", "SD02", "SD03")]
#data filtering------
# filtering with out all genes max count < 50.
# any gene with maximum count lower than 50 will be filtered out
keep <- apply(counts, 1, max) >= 50
counts <- counts[keep,]
dim(counts)  # save this number in your notebook. as it is the number of genes that will be used for the downstream analysis


# Defining a design matrix based on the experimental design----
design <- matrix(data=c(1,1,1,0,0,0,0,0,0,1,1,1), nrow=6, ncol=2, dimnames = list(c(), c('DO','SD')))
design
#creating a matrix of contrasts, where which each column represents a contrast between two groups of interest.
cont.matrix <- matrix(data=c(-1,1), nrow=2, ncol=1, dimnames = list(c(), c('SD')))
cont.matrix
# Create a DGEList object-----
y <- DGEList(counts=counts)
#normalization----
#factor normalization
y <- calcNormFactors(y, method="TMM")
#estimate dispersion-----
y <- estimateDisp(y, design, robust=TRUE)

#statistics tests------
# Now conduct glm fit test for treatment effect
fit <- glmFit(y,design)
lrt <- glmLRT(fit, contrast=cont.matrix)
#The total number of genes significantly up-regulated or down-regulated at 5% FDR is summarized as follows:
summary(decideTests(lrt, p.value = 0.01, lfc=1))

# the top 10 DEGs----
topTags(lrt)
FDR <- p.adjust(lrt$table$PValue, method="BH")
sum(FDR < 0.05)

#MA (mean average plot) plot-----
plotMD(lrt, main="IMtS vs RX")
abline(h=c(-1,1), col="blue")
#generate the table of results and save it in a file
out <- topTags(lrt, n=Inf, sort.by='PValue')$table

#merging------
#merge the statistics results with the feature_count and the annotation table
m1 <- merge(out, counts, by.x=0,by.y=0)
m2 <- merge(m1,gp, by.x='Row.names', by.y='Geneid')
dim(m2)
dim(m1)
dim(out)
#save it it in a file 
fwrite(m2, file="stat1_SD_vs_DO.csv")

#volcano_plots
# how improve this visualization
m2$FDR[m2$FDR <1e-20] <- 1e-40
ggplot(data=m2, aes(x=logFC, y=-log10(FDR), col=Target))+geom_point(stroke =2 , shape = 16, size=3)+
  theme(panel.background = element_rect(fill = "white", color = "grey50"))+
  scale_color_manual(breaks = c("COEs", "CUTs", "CYPs", "GSTs","SGPs","OBPs", "Others"),
                     values=c("red", "green", "blue","violet", "black","pink","darkgrey"))+
  #scale_y_continuous(limits = c(0, 120))+
  geom_vline(xintercept=c(-1, 1), col="black", lty=3) + geom_hline(yintercept = -log10(0.01), lty=3, color="black")+
  theme(legend.title=element_blank())+
  #scale_y_continuous(limits = c(0, 90))+
  ggtitle("I-MAL vs Rock")+
  labs(x = expression ("log"[2](FC)), y=expression("-log"[10](FDR)))