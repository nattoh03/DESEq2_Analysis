# DESEq2_Analysis
Do different sample groups express genes differentially



Three questions to answer


#### which genes do the reads (samples) belong to ?
    done in supercomputer, alignment of ref genome with Hisat2

### How many reads align to a specifi gene?
   done in supercomputer, using htseq-count, the gff file and the sorted.bam file to generate a txt/csv file 
   
##### Do different samples express genes differentially?
    now this part is done in R using library(DESeq2)
    
    

    
