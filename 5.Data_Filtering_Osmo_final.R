#################################################################
### load count data or depth data (with or without row.names) ###
#################################################################

# load count data
countData <- read.csv("genecount.csv", header=TRUE, row.names=1)

# check the number of samples/genes
dim(countData)

### remove rows with 0 counts  ###
countData$missing_number <- rowSums(countData == 0)

sum(countData$missing_number) #check if there are any that are not expressed in all samples

# histogram with missing gene expression

pdf("geneCounts_missing_hist.pdf")
hist(countData$missing_number, col="deepskyblue", border=F)
dev.off()

# check the missing genes
summary(countData$missing_number)

countData_missing_removed <- subset(countData, missing_number <= 3) # genes expressed by at least 50% of the samples; relax the threshold if needed

length(countData$missing_number)-length(countData_missing_removed$missing_number)  

# check how many genes were removed
dim(countData_missing_removed)

# delete the column with the missingness
countData_missing_removed$missing_number <- NULL

# change the counts to numeric; maybe this needs to be changed to integers
all(is.numeric(countData_missing_removed))
#countData_missing_removed <-apply(countData_missing_removed, 2, as.numeric)

write.csv(countData_missing_removed, file="GeneCount_Clean_Gene_Names.csv", row.names = TRUE)

