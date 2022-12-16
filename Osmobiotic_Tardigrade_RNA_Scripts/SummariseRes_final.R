###################################################
# Summarise the results from the three DE methods #
###################################################

# R
# remove anything before we start
rm(list = ls())

# Read files
# Results files from the three methods Limma, EdgeR and DESeq2:
m1<-read.csv("Limma_res_Ctl24h_Exp24h_final.csv")
m2<-read.csv("EdgeR_glmQLF_res_Ctl24h_Exp24h.csv")
m3<-read.csv("DESeq2_res_Ctl24h_Exp24h_Osmobiosis_final.csv")

# read the annotation files
annotation_info1<-read.table("Ro-001-Unigene.fa.blast.swissprot.txt", as.is=T, sep="\t", h=T)
annotation_info2<-read.table("Ro001UnigeneGene2GOedited", as.is=T, sep="\t")

# the 1 inside the brackets correspond to the column that contain the gene id (unigene or contig ID), so if it is not column 1 then change accordingly
DE_m1<-unique(m1[,1],fromLast = TRUE, nmax = NA)
DE_m2<-unique(m2[,1],fromLast = TRUE, nmax = NA)
DE_m3<-unique(m3[,1],fromLast = TRUE, nmax = NA)

DE_genes<-unique(c(DE_m1,DE_m2,DE_m3))

m1_p<-NULL
m2_p<-NULL
m3_p<-NULL
gene_id<-NULL
GO_id<-NULL
descrip<-NULL
m1_FDR<-NULL
m2_FDR<-NULL
m3_FDR<-NULL

for(i in 1:length(DE_genes)){
  # change here the 1 if that column doesn't corresponds to the gene id
  if(sum(m1[,1]==DE_genes[i])>0){
    # Here you could change "method 1" for the name of the program you used (e.g. limma, DESeq)
    m1_p<-c(m1_p, "Limma_Ctl24h_Exp24h")
    # Here the 1 in "m1[,1]" should be the column that contains the gene id, and the 2 in "[i],2]" should be the column that contains the FDR
    m1_FDR<-c(m1_FDR, m1[m1[,1]==DE_genes[i],7])
  }else{
    m1_p<-c(m1_p, NA)
    m1_FDR<-c(m1_FDR, NA)
  }
  
  # change here the 1 if that column doesn't corresponds to the gene id
  if(sum(m2[,1]==DE_genes[i])>0){
    # Here you could change "method 2" for the name of the program too
    m2_p<-c(m2_p, "EdgeR_Ctl24h_Exp24h")
    # Here the 1 in "m2[,1]" should be the column that contains the gene id, and the 2 in "[i],2]" should be the column that contains the FDR
    m2_FDR<-c(m2_FDR, m2[m2[,1]==DE_genes[i],6])
    
  }else{
    m2_p<-c(m2_p, NA)
    m2_FDR<-c(m2_FDR, NA)
  }
  
  # change here the 1 if that column doesn't corresponds to the gene id
  if(sum(m3[,1]==DE_genes[i])>0){
    # Here you could change "method 1" for the name of the program too
    m3_p<-c(m3_p, "DESeq2_Ctl24h_Exp24h")
    # Here the 1 in "m3[,1]" should be the column that contains the gene id, and the 2 in "[i],2]" should be the column that contains the FDR
    m3_FDR<-c(m3_FDR, m3[m3[,1]==DE_genes[i],7])
    
  }else{
    m3_p<-c(m3_p, NA)
    m3_FDR<-c(m3_FDR, NA)
  }
  
  if(sum(annotation_info1[,1]==DE_genes[i])>0){
    temp<-annotation_info1[annotation_info1[,1]==DE_genes[i],2]
    gene_id<-c(gene_id, temp[1])
    temp<-annotation_info1[annotation_info1[,1]==DE_genes[i],13]
    descrip<-c(descrip, temp[1])
  }else{
    gene_id<-c(gene_id, NA)
    descrip<-c(descrip, NA)
  }
  if(sum(annotation_info2[,1]==DE_genes[i])>0){
    GO_id<-c(GO_id, annotation_info2[annotation_info2[,1]==DE_genes[i],2])
  }else{
    GO_id<-c(GO_id, NA)
  }	
}

tab<-cbind(DE_genes, m1_p, m1_FDR, m2_p, m2_FDR, m3_p, m3_FDR, gene_id, GO_id, descrip)

colnames(tab)<-c("Seq_id", "m1", "m1_FDR", "m2", "m2_FDR", "m3", "m3_FDR", "gene_id", "GO_kegg", "Description")

#It gets ordered strangely 
write.csv(tab, "Summary_table_csv.csv")

Sum <- read.csv("Summary_table_csv.csv")
Sum_GeneID <- Sum[,c(2,9)]
na_vec <- which(!complete.cases(Sum_GeneID))

ListOfGenes <- Sum_GeneID[-na_vec,]

ListOfGeneID <- ListOfGenes[2]

write.csv(as.data.frame(ListOfGeneID), file="ListOfGeneID.csv")

####################################



