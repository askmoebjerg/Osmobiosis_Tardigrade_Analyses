###############################################
# Enrichment analysis using topGO library

# remove anything before we start
rm(list = ls())

library(GO.db)
library(topGO)
library(Rgraphviz)

# set the working directory
setwd("C:/Users/askmo/OneDrive - Københavns Universitet/Uni MolBio/Articles I coauthor/Osmotic Stress in Tardigrades")

# read Ro001UnigeneGene2GOedited which contains the annotation data 
geneID2GO <- readMappings("Ro001UnigeneGene2GOedited", sep = "\t", IDsep = ",")
str(head(geneID2GO))

# Define a list of DE genes #
geneNames <- names(geneID2GO)
head(geneNames)

m1<-read.csv("Limma_res_Ctl24h_Exp24h_final.csv")
m2<-read.csv("EdgeR_glmQLF_res_Ctl24h_Exp24h.csv")
m3<-read.csv("DESeq2_res_Ctl24h_Exp24h_Osmobiosis_final.csv")

d <- read.csv("Summary_table_csv.csv")
n123 = subset(d, m1=="Limma_Ctl24h_Exp24h" & m2=="EdgeR_Ctl24h_Exp24h" & m3=="DESeq2_Ctl24h_Exp24h", select = Seq_id)
DE_across <- unique(n123[,1],fromLast = TRUE, nmax = NA)

#For investigating DEGs
GenesOfInterest <- DE_across

geneList <- factor(as.integer(geneNames %in% GenesOfInterest)) #binary list indicating the DE genes
names(geneList) <- geneNames
str(geneList)

# create an topGO data object #
myGOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

# run the Fisher's exact tests
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
resultElim <- runTest(myGOdata, algorithm="elim", statistic="fisher")
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
resultParentchild <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")

# see how many results we get where weight01 gives a P-value <= 0.05:
mysummary <- summary(attributes(resultTopgo)$score <= 0.05)
numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.05

# create and write a table summarising the top 'numsignif' results:
allRes <- GenTable(myGOdata, classicFisher = resultClassic, elimFisher = resultElim, topgoFisher = resultTopgo, parentchildFisher = resultParentchild, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif)

write.csv(allRes, "TopGO_Results_DE_across.csv")

# print a graph (to a pdf file) with the top 'numsignif' results:
printGraph(myGOdata, resultTopgo, firstSigNodes = numsignif, fn.prefix = "topGO", useInfo = "all", pdfSW = TRUE)
dev.off()

