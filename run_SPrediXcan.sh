#!/bin/bash

#PAGE data
#use Ryan's models, baseline
#/home/isabelle/topmed/proteome/PAGE

p="ALL"
#p2="PAV_filtered_rho0.1_zpval0.05"
#p2="PAV_filtered_unfiltered"
#p2="rho0.1_zpval0.05"
p2="unfiltered"

FILES=/home/elyse/sumstats_for_topmed/Wojcik_build38/WojcikG_*.txt.gz
for f in $FILES
do
  /usr/local/bin/MetaXcan_software/SPrediXcan.py \
  --model_db_path /home/ryan/TOPMed_Proteome/dbs_out/${p}_PCAIR_baseline_models_${p2}.db \
  --covariance /home/ryan/TOPMed_Proteome/dbs_out/${p}_PCAIR_baseline_models_${p2}_covariances.txt.gz \
  --gwas_file $f \
  --snp_column SNP_hg38 \
  --effect_allele_column Effect-allele \
  --non_effect_allele_column Other-allele \
  --beta_column Beta \
  --se_column SE \
  --pvalue_column P-val \
  --output_file /home/isabelle/topmed/proteome/PAGE/${p}_PCAIRbaseline_${p2}_${f:55:(-7)}.csv \
  --keep_non_rsid
done


#string manipulation in BASH https://www.tldp.org/LDP/abs/html/string-manipulation.html
#used to rm .txt.gz from output file name
#also the path from the file name
