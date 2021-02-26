#this script takes the [model]_bonferroni_all_pheno.csv file and pulls unique gene-trait pairs, generates list of JUST gene_name and Phenotype [model]_unique_sig_all_pheno.csv
library(dplyr)
library(data.table)
"%&%" = function(a,b) paste(a,b,sep="")

model_types<-c('PAV_filtered_rho0.1_zpval0.05','PAV_filtered_unfiltered','rho0.1_zpval0.05','unfiltered')
for (model in model_types){
	data<-read.csv('PCAIRbaseline_'%&%model%&%'_bonferroni_all_pheno.csv')
	unique_genepairs <- data %>% select(gene_name, Phenotype)
        unique_genepairs <- unique(unique_genepairs)
        print(c(model, 'unique significant gene trait pairs: ', nrow(unique_genepairs)))
	write.csv(unique_genepairs,'PCAIRbaseline_'%&%model%&%'_unique_sig_all_pheno.csv',quote=F,row.names=F)
}
