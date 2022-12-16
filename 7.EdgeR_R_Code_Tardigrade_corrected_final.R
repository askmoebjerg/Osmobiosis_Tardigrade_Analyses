#5.EdgeR_R_Code_Tardigrade_4MAr22.R

# remove anything before we start
rm(list = ls())

# Install packages
library(DESeq2)
library(data.table)
library(edgeR)
library(limma)
library(xlsx)
library(statmod)
library(DEFormats)

############################
### load phenotypic data ###
############################

# Load phenotypic data: should contain sample name, phenotype, date of extraction, person who extracted the RNA, lane used for sequencing...etc
colData <- read.csv("Samplename_Osmotolerance.csv", header=TRUE)

# Format phenotypic data
colData$SampleID <- as.factor(colData$SampleID)
colData$Group <- as.factor(colData$Group) # control or experiment
colData$Thawing_Date <- as.factor(colData$Thawing_Date)
colData$Extraction_Date <- as.factor(colData$Extraction_Date)
colData$Extractor <- as.factor(colData$Extractor)
colData$Nr_Specimens <- as.integer(colData$Nr_Specimens)

#################################################################
### load count data or depth data (with or without row.names) ###
#################################################################

# Load a count matrix
countData <- read.csv("GeneCount_Clean_Gene_Names.csv", header=TRUE, row.names=1)


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

###################
### Start EdgeR ###
###################

### on the whole countData ###
x <- DGEList(counts=countData, group=colData$group)

# Normalize for RNA composition
x <- calcNormFactors(x)
x$samples

# Compute counts per million (CPM)
cpm <- cpm(x)

# MDS (Plot samples on a two-dimensional scatterplot)
#plotMDS(x, col=group)
plotMDS(x,cex=0.5)

# Design matrix
design <- model.matrix(~ 0 + group, data=colData)
design
colnames(design) # to get the group names to use for contrasts

# Estimate the dispersions
x <- estimateDisp(x, design, robust=TRUE)

#Scatterplot of the biological coefficient of variation (BCV) against the average abundance of each gene.
#stabilisations of the variables 
plotBCV(x)


################################################
##### To perform quasi-likelihood F-tests ######
################################################

fit <- glmQLFit(x, design, robust=TRUE)
plotQLDisp(fit)

# Smaller prior df estimates indicate that the true unknown dispersions are highly variable, so weaker moderation towards the trend is appropriate.
summary(fit$df.prior)

##### Make contrasts #####

# Ctl_24h vs. Exp_24h - 3372 Down, 2459 UP, 5831 Total DE genes
Ctl24h_Exp24h <- makeContrasts(groupExp_Osmobiosis_24h-groupCtl_Osmobiosis_24h, levels=design)
qlf_Ctl24h_Exp24h  <- glmQLFTest(fit, contrast=Ctl24h_Exp24h)
is.de <- decideTestsDGE(qlf_Ctl24h_Exp24h, adjust.method="BH", p.value=0.05)
summary(is.de)
out_Ctl24h_Exp24h <- topTags(qlf_Ctl24h_Exp24h, n=Inf, adjust.method="BH")
keep_Ctl24h_Exp24h <- out_Ctl24h_Exp24h$table$FDR <= 0.05
output_Ctl24h_Exp24h <- out_Ctl24h_Exp24h[keep_Ctl24h_Exp24h,]
write.csv(as.data.frame(output_Ctl24h_Exp24h), file="EdgeR_glmQLF_res_Ctl24h_Exp24h_corrected.csv")


#########################################################################
##### Creating log-fold change against log-counts per million plot ######
#########################################################################

#Plot log-fold change against log-counts per million, with DE genes highlighted:
dev.off()
plotMD(qlf_Ctl24h_Exp24h, p.value = 0.05, main = "MD-Plot", cex=0.35)
legend("topright",c("NotSig","Up","Down"),col=c("Black","red","blue"),pch=16, cex = 0.9)



#############################################
##### To perform likelihood ratio tests #####
#############################################

fit_2 <- glmFit(x, design, robust=TRUE)
lrt_2 <- glmLRT(fit_2)

# Ctl_24h vs. Exp_24h - 2733 Down, 2192 UP, 4926 Total  DE genes
Ctl24h_Exp24h <- makeContrasts(groupExp_Osmobiosis_24h-groupCtl_Osmobiosis_24h, levels=design)
lrt_Ctl24h_Exp24h  <- glmLRT(fit, contrast=Ctl24h_Exp24h)
is.de <- decideTestsDGE(lrt_Ctl24h_Exp24h, adjust.method="BH", p.value=0.05)
summary(is.de)
out_Ctl24h_Exp24h <- topTags(lrt_Ctl24h_Exp24h, n=Inf, adjust.method="BH")
keep_Ctl24h_Exp24h <- out_Ctl24h_Exp24h$table$FDR <= 0.05
output_Ctl24h_Exp24h <- out_Ctl24h_Exp24h[keep_Ctl24h_Exp24h,]
write.csv(as.data.frame(output_Ctl24h_Exp24h), file="EdgeR_glmLRT_res_Ctl24h_Exp24h_corrected.csv")


#########
# Limma #
#########

dge <- DGEList(counts=countData)

design <- model.matrix(~ 0 + group, data=colData)

dge <- calcNormFactors(dge)

v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)
head(coef(fit))


### Ctl_24h vs. Exp_24h - 4446 Down, 3290 Up, 7737 Total  DE genes ###
contr <- makeContrasts(groupExp_Osmobiosis_24h - groupCtl_Osmobiosis_24h, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
summary(decideTests(tmp))
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- top.table[ which(top.table$adj.P.Val < 0.05), ]
write.csv(top.table, file = "Limma_res_Ctl24h_Exp24h_corrected.csv")

summary(decideTests(tmp))

