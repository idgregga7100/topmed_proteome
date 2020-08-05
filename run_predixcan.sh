#!/bin/bash
#PBS -S /bin/bash
#PBS -o /home/isabelle/topmed/proteome/logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e /home/isabelle/topmed/proteome/logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
#PBS -l walltime=150:00:00
#PBS -l mem=16gb

p="ALL"
p2="0.01"
p3="T"

/usr/local/bin/MetaXcan_software/PrediXcan.py \
--throw \
--text_genotypes /home/isabelle/topmed/proteome/geno_dosages/no-header_chr*.maf0.01.dosage.txt.gz \
--text_sample_ids /home/isabelle/topmed/proteome/sample_list.txt \
--model_db_path /home/ryan/topmed/proteome/dapg_net_CORRECT_WINDOWS/current_dbs/${p}_separate_WG_PIP_${p2}_cluster_filt_${p3}_unfiltered.db \
--prediction_output /home/isabelle/topmed/proteome/results/${p}_PIP_${p2}_${p3}_model_predicted_expression.txt \
--prediction_summary_output /home/isabelle/topmed/proteome/results/${p}_PIP_${p2}_${p3}_model_prediction_summary.txt

echo /home/isabelle/topmed/proteome/results/${p}_PIP_${p2}_${p3}_model_prediction_summary.txt
