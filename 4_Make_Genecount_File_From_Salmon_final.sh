#############################
# Creating gene count file  #
#############################

cut -f1 T1_quant.sf > genename
for f in T1_quant.sf T2_quant.sf T3_quant.sf T4_quant.sf T5_quant.sf T6_quant.sf
do
    sample=$f
    awk 'BEGIN {OFS="\t"} {print $5}' $f | sed -e "1s/NumReads/$sample/" > $f.tmp
done
paste genename T1_quant.sf.tmp T2_quant.sf.tmp T3_quant.sf.tmp T4_quant.sf.tmp T5_quant.sf.tmp T6_quant.sf.tmp > genecount

########################################
# convert the gene count table to .csv #
########################################

cat genecount | tr -s '[:blank:]' ',' > genecount.csv

# remove trailing "quant.sf" text  
sed 's/_quant.sf//' genecount.csv > genecount.csv
