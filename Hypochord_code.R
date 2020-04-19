library(edgeR)
library(DESeq2)
library(limma)
library(Biobase)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library("org.Mm.eg.db")
library(RColorBrewer)
library(gplots)
library(knitr)
library(calibrate)
library(pheatmap)
library(RColorBrewer)
library(gProfileR)
library("ggplot2")
library("ggdendro")
library("reshape2")
library("grid")
library("limma")
library("MeSH.Xtr.eg.db")
library("dplyr")
library("ggplot2")
library("PoiClaClu")


seqdata <- read.table("count_matrixH.tabular") 
#removing first columns from the uploaded file
countdata <- seqdata[,-(1)] 
rownames(countdata) <- seqdata[,1]
head(countdata)

#defining the conditions for the hypochord
condition <-factor(c (rep("methimazole", 1), rep("methimazole", 1), rep("control", 1), 
                      rep("control", 1), rep("control", 1)))
timepoints <-factor(c (rep ("t1", 2), rep("t3", 1), rep("t2", 1), rep("t1", 1)))
sampleTable <-data.frame(timepoints = as.factor(timepoints), condition = as.factor(condition))
rownames(sampleTable)=colnames(countdata)
sampleTable

#PCA plot for the log transformed data
rld <- rlog(dds, blind = FALSE)
plotPCA(rld, intgroup = c("condition", "timepoints"))

#Differential expression analysis using DESeq2

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = sampleTable, design = ~timepoints+condition)
dds

#remove counts with less than 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
head(dds)
dds$condition <- factor(dds$condition, levels = c("methimazole","control"))
dds$timepoints <-factor(dds$timepoints, levels = c("t1","t2", "t3"))




#running dds
dds <- DESeq(dds)


#Differential expression analysis
#Specifying the coefficient or contrast we want to build a results table for...
#Following command gives the contrast states you can test for the differential expression analysis
resultsNames(dds) 
res <- results(dds, contrast=c("condition","methimazole","control"))
resultsNames(dds) #this gives the contrasts that can be done
res
summary(res)
dds <- DESeq(dds, betaPrior = FALSE)



