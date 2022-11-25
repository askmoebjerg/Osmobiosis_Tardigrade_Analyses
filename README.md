# Osmobiosis_Tardigrade_Analyses
The following contains a collection of scripts used to perform differential gene expression analyses on Tardigrade RNA data derived from osmobiosis experiments. 

## Pipeline
The scripts used for the analyses can be found in the `Osmobiotic_Tardigrade_RNA_Scripts` 
- Initial quality assessment of the sesuencing data by running FastQC on the "CLEAN" reads
- Run Salmon to quantify the number of read counts 
- Creating a ´CSV´ input file from the individual gene count file
- Filtering out the genes with too many missing counts 
- Perform differential gene expression (DE) analyses 
      -  DE using DESeq2
      -  DE using EdgeR and Limma 
- Creating a summary table with the differentially expressed genes found by DESeq2, EdgeR and Limma 
- Perform the Gene Ontology (GO) enrichment analyses on the identified differentially expressed genes    
