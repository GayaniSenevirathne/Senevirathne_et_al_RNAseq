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
library("EnhancedVolcano")
library("tidyverse")
library("dplyr")
library("ggplot2")


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
res1 <- results(dds, contrast = c("condition","methimazole","control"))
res1 <- lfcShrink(dds, contrast = c("condition","methimazole","control"), res = res1)

resOrdered <- res1[order(res1$pvalue),]
summary(resOrdered)
summary(res1)

#gives the numbers of genes which are at pvalue o.1 or 0.05
sum(res1$padj < 0.1, na.rm=TRUE)
sum(res1$padj < 0.05, na.rm=TRUE)

#alpha defines the adjusted pvalue == FDR value)
res05 <- results(dds, contrast=c("condition","methimazole","control"), alpha=0.05)
summary(res05)

#subset the results table to these genes and then sort it by the log2 fold change estimate to get the significant genes,
#with the strongest down-regulation:
resSig <- subset(res, padj < 0.05) #lower the pvalue, better. usually 0.05 can be taken.
summary(resSig)

#getting gene annotations from biomart
#this could be done for any dataset: either the "res" or "vsd" BUT..
#here we are annotating resSig,  which is a list of genes taken from "urostylevshypochord"
library(biomaRt)
ensembl <- useMart('ensembl', 
                   dataset = 'xtropicalis_gene_ensembl',
                   host = "www.ensembl.org")

filterValues <- rownames(resSig)
genes <- getBM(mart = ensembl, attributes = c('ensembl_gene_id',
                                              'external_gene_name'), filters = c("ensembl_gene_id"), values = filterValues)

#remove duplicates: 
genes <- genes[!duplicated(genes$ensembl_gene_id),]


#annotating resSig from the biomart dataset
res_annot <- as.data.frame(resSig) %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  left_join(genes, "ensembl_gene_id") 
res_annot
all(rownames(resSig) == res_annot$ensembl_gene_id)


#getting up and down regulatory genes, between urostyle and hypochord.

down_new<-res_annot[ order(resSig$log2FoldChange,-resSig$baseMean), ]
write.csv(as.data.frame(down_new), 
          file="results_alpha_0.05_hypochord_Meth_vs_control_downregulating.csv")


#and the genes with strongest upregulation
Up_new<-res_annot[ order(-resSig$log2FoldChange,
                         -resSig$baseMean), ]
write.csv(as.data.frame(Up_new), 
          file="results_alpha_0.05_hypochord_Meth_vs_control_upregulated_genes.csv")


##heatmap for selected genes
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

mat  <- assay(vsd)
head(rownames(mat))
idx <- rownames(mat)[ which(res$padj < 0.05)]
mat[ idx, ]


filterValues <- rownames(mat[ idx, ])
genes <- getBM(mart = ensembl, attributes = c('ensembl_gene_id',
                                              'external_gene_name'), filters = c("ensembl_gene_id"), values = filterValues)

head(genes)

#remove duplicates: 
genes <- genes[!duplicated(genes$ensembl_gene_id),]

#annotating resSig from the biomart dataset
mat_anno <- as.data.frame(mat[ idx, ]) %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  left_join(genes, "ensembl_gene_id") 
mat_anno
all(rownames(mat[ idx, ]) == mat_anno$ensembl_gene_id)

#drawing the heatmap for the selected genes from the DE genes

mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c("purple","orange")[sampleTable$condition] #color given according to condition (urostyle/hypochord)
col.status <- c("blue","red","dark green")[sampleTable$timepoints] #color given according to the timepoints
library(devtools)
library(ggplot2)
clab<-cbind(col.cell,col.status) #bind the two color formats together

#get the heatmap.3 function downloaded
library(devtools)
library(ggplot2)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

pdf(file="p_0.0005_hypochord_vs_urostyle_geneIDs_50_genes.pdf")
heatmap.3(mat[ idx, ],col=rev(morecols(50)),trace="none", main="Hypochord_DEGs_MethimazoleVSControl", ColSideColors=clab, scale="row", labRow=mat_anno$external_gene_name)
dev.off()




#volcanoplot for selected genes

keyvals <- rep("black", nrow(res1))
names(keyvals) <- rep("Mid / NA", nrow(res1))

keyvals[which(res1$log2FoldChange > 7)] <- "purple"
names(keyvals)[which(res1$log2FoldChange > 7)] <- "High"

keyvals[which(res1$log2FoldChange < -7)] <- "gold"
names(keyvals)[which(res1$log2FoldChange < -7)] <- "Low"


pdf(file="Hypo_vs_Meth_volcanoplot_Selected_genes_labelled2.pdf")  

EnhancedVolcano(res1,
                lab = rownames(res1),
                x = "log2FoldChange",
                y = "padj",
                selectLab = rownames(res1)[which(names(keyvals) %in% c("High", "Low"))],
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1.5,
                xlim = c(-15,15),
                transcriptLabSize = 2.0,
                colAlpha = 1,
                legendPosition = "top",
                transcriptPointSize = 1.8,
                legendLabSize = 6,
                legendIconSize = 4.0,
                DrawConnectors = TRUE,
                widthConnectors = 0.3,
                colConnectors = "grey50",
                border = "partial",
                borderWidth = 1.5,
                borderColour = "black",
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                colOverride = keyvals)
dev.off()
