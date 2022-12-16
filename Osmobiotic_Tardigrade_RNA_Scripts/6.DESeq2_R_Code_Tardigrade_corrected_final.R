#4.DESeq2_R_Code_Tardigrade_corrected.R

# remove anything before we start
rm(list = ls())

library(DESeq2)
library(data.table)
library(edgeR)
library(limma)
library(xlsx)
library(statmod)
library(DEFormats)
library(sva)

############################
### load phenotypic data ###
############################

# Load phenotypic data
colData <- read.csv("Samplename_Osmotolerance.csv", header=TRUE)

# Format phenotypic data
colData$SampleID <- as.factor(colData$SampleID)
colData$Group <- as.factor(colData$Group) # control or experiment
colData$Thawing_Date <- as.factor(colData$Thawing_Date)
colData$Extraction_Date <- as.factor(colData$Extraction_Date)
colData$Extractor <- as.factor(colData$Extractor)
colData$Nr_Specimens <- as.integer(colData$Nr_Specimens)

#######################
### load count data ###
#######################

# Load a count matrix
countData <- read.csv("GeneCount_Clean_Gene_Names.csv", header=TRUE, row.names=1)

####################
### Start DESeq2 ###
####################

# Change the count table to integer
countData[] <- lapply(countData, as.integer)

#correct for batch effect 
countData <- sva::ComBat_seq(as.matrix(countData), batch = colData$Thawing_Date,shrink=FALSE, shrink.disp=FALSE)

# Full data - make groups, and put that into "design"
colData$group <- factor(paste0(colData$Group, "_", colData$Duration))

#For Control comparison
colData$group <- factor(paste0(colData$Group, "_", colData$Experiment, "_", colData$Duration))

#this command allows linear regression, we now choose the group with small g
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ group)
model.matrix(~ colData$group)
dds$group<-relevel(dds$group,ref="Ctl_Osmobiosis_24h") #To set the control as reference

### 4.2The variance stabilizing transformation and the rlog ###
# The VST is much faster to compute and is less sensitive to high count outliers than the rlog. The rlog tends to work well on small datasets (n < 30), potentially outperforming the VST when there is a wide range of sequencing depth across samples (an order of magnitude difference). We therefore recommend the VST for medium-to-large datasets (n > 30).
#removes outliers 

vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)

library("dplyr")
library("ggplot2")

### 4.3 Sample distances ###

# calculate the Euclidean distance between samples. To ensure we have a roughly equal contribution from all genes, we use it on the VST data.
sampleDists <- dist(t(assay(vsd)))
sampleDists

library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(vsd$Group, vsd$Experiment, vsd$Duration, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

### 4.4PCA plot
#this allows to us to see variations 
# plots can be customized using ggplot function too
plotPCA(vsd, intgroup=c("Group", "Experiment", "Duration"))
plotPCA(vsd, intgroup = "Thawing_Date")
plotPCA(vsd, intgroup = "Extraction_Date")
plotPCA(vsd, intgroup = "Extractor")

### 5.1Running the differential expression pipeline
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
summary(res)

###

# check levels
levels(dds$group)

######

# gene with lowest p-value
plotCounts(dds, gene=which.min(res$padj), intgroup="group")

#######

# Ctl_24h vs. Exp_24h - total 3346, down 1687, up 1659 DE genes
res_Ctl24h_Exp24h <- results(dds, contrast=c("group", "Exp_Osmobiosis_24h","Ctl_Osmobiosis_24h"), alpha=0.05) 
summary(res_Ctl24h_Exp24h)
sum(res_Ctl24h_Exp24h$padj < 0.05, na.rm=TRUE)
res_Ctl24h_Exp24h <- res_Ctl24h_Exp24h[ which(res_Ctl24h_Exp24h$padj < 0.05), ]
res_Ctl24h_Exp24h_ordered <- res_Ctl24h_Exp24h[order(res_Ctl24h_Exp24h$padj),]

write.csv(as.data.frame(res_Ctl24h_Exp24h_ordered), file="DESeq2_res_Ctl24h_Exp24h_Osmobiosis_corrected.csv") 

plotCounts(dds, gene=which.min(res_Ctl24h_Exp24h$padj), intgroup="group")

#########################################################################
##### Creating log-fold change against log-counts per million plot ######
#########################################################################

#Plot log-fold change against log-counts per million, with DE genes highlighted:
dev.off()
order_res_mean <- res_Ctl24h_Exp24h[order(res_Ctl24h_Exp24h$baseMean),]
 
DESeq2::plotMA(order_res_mean, main = "MA-Plot Osmobiosis 24h", ylim = c(-4, 4), hl.col= c("blue","red"))

resG <- results(dds, lfcThreshold=0, altHypothesis="greater")
resL <- results(dds, lfcThreshold=0, altHypothesis="less")
resG <- results(dds, contrast=c("group", "Exp_Osmobiosis_24h","Ctl_Osmobiosis_24h"), alpha=0.05,altHypothesis="greater") 

DESeq2::plotMA(resG, ylim=c(-4, 4),)
points(resL$log2FoldChange)
DESeq2::plotMA(resL, ylim=c(-4, 4), sigCol)


############################# 
# Check Ctl_24h vs. Exp_24h #
# DE genes - 3302           #
#############################

# Change the count table to integer
countData[] <- lapply(countData, as.integer)

# Full data - make groups, and put that into "design"
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Group)
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)

