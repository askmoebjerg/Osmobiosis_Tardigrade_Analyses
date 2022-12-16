##############
# Run Salmon #
##############

#############################################
# Index see DOI: 10.1016/j.cbpa.2022.111245 #
#############################################

# Run salmon
# -l A tells salmon that it should automatically determine the library type of the sequencing reads (e.g. stranded vs. unstranded etc.).
for SAMPLE in T1 T2  T3  T4  T5  T6
do
    salmon quant -i /corefac/cbp/tardigrades_project/analysis/Salmon/Salmon_tardigrade_bgi_de_novo_index -l A -1 /corefac/cbp/tardigrades_project/BGI_data/$SAMPLE/*_1.fq -2 /corefac/cbp/tardigrades_project/BGI_data/$SAMPLE/*_2.fq -p 2 -o $SAMPLE"_quant"
done
