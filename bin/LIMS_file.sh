#!/usr/bin/env bash
# This script takes the outputs form abricate and merges them into one tabsepearted file for imprting into LIMS casebook

# Extract data from *_sero.csv files
awk 'FNR==1 && NR!=1 { while (/^#F/) getline; } 1 {print}' *_serotype.csv > sero_file.csv

# Add "CATEGORY" header to sero_file.csv
awk 'BEGIN{print "CATEGORY"} {if(NR>1) print ($0=="" ? "" : "serotype")}' sero_file.csv > sero_column.txt

# Combine sero_file.csv and sero_column.txt into serotype_res.csv
paste sero_file.csv sero_column.txt > serotype_res.csv

# Extract data from *_sero.csv files
awk 'FNR==1 && NR!=1 { while (/^#F/) getline; } 1 {print}' *sero_Jia.csv > jia_sero_file.csv


# Extract data from *_vf.csv files
awk 'FNR==1 && NR!=1 { while (/^#F/) getline; } 1 {print}' *_vf.csv > vf_file.csv

# Add "CATEGORY" header to vf_file.csv
awk 'BEGIN{print "CATEGORY"} {if(NR>1) print ($0=="" ? "" : "VF")}' vf_file.csv > vf_column.txt

# Combine vf_file.csv and vf_column.txt into vf_res.csv
paste vf_file.csv vf_column.txt > vf_res.csv

# Extract data from *_AMR.csv files
awk 'FNR==1 && NR!=1 { while (/^#F/) getline; } 1 {print}' *_AMR.csv > amr_file.csv

# Add "CATEGORY" header to amr_file.csv
awk 'BEGIN{print "CATEGORY"} {if(NR>1) print ($0=="" ? "" : "AMR")}' amr_file.csv > amr_column.txt

# Combine vf_file.csv and vf_column.txt into amr_res.csv
paste amr_file.csv amr_column.txt > amr_res.csv

# Generate current date and time
datetime=$(date +"%d%b%Y_%H-%M-%S")



# Extract data from *_res.csv files and save as ${datetime}_LIMS_file.csv
awk 'FNR==1 && NR!=1 { while (/^#F/) getline; } 1 {print}' *_res.csv > draft_file.tsv

cat software_version.tsv draft_file.tsv >> ${datetime}_Gpara_LIMS_file.csv

# Replace "#FILE" with "ID" in ${datetime}_LIMS_file.csv
sed -i 's,#FILE,ID,g' ${datetime}_Gpara_LIMS_file.csv
sed -i 's,_assembly.fasta,,g' ${datetime}_Gpara_LIMS_file.csv


