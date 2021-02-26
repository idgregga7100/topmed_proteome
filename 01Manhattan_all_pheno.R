#elyse's script /home/egeoffroy/topmed/SPrediXcan/scripts/Manhattan_all_pheno.R
#had to edit it to work for my directories
#this entire script generates two output files: a [model]_sig_genes_all_pheno.csv file and a [model]Manhattan_total_phenotype_5e8.tiff file
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(qqman))
suppressPackageStartupMessages(library(wesanderson))
options(warn=-1)

#need to make BP Chrome files to use here, used elyse's make_BP_Chrome_2.py
load_bp_chrom <- function(pop, type){
        print(paste("/home/isabelle/topmed/proteome/PAGE/BP_Chrom_files/", pop, "_", type, "_BP_Chrome.csv", sep = ''))
        BP_Chrome <- read.table(paste("/home/isabelle/topmed/proteome/PAGE/BP_Chrom_files/", pop, '_', type, "_BP_Chrome.csv", sep = ''), sep = ',', header = T)
        colnames(BP_Chrome) <- c('PROBE_ID', 'gene_name', 'chr', 'BP')

        BP_Chrome <- BP_Chrome %>%  
                transform(chr = str_replace(chr, "chr", "")) 
        BP_Chrome <- transform(BP_Chrome, chr=as.numeric(chr))
        BP_Chrome <- transform(BP_Chrome, BP=as.numeric(BP))
        BP_Chrome$GENE <- gsub("\\..*","",BP_Chrome$PROBE_ID)
        return(BP_Chrome)
}

get_files <- function(model_type){
	files <- list.files('/home/isabelle/topmed/proteome/PAGE/SPred_out', pattern = model_type, recursive = T, full.names=T)
	return(files)
}

format_files <- function(files, model_type, output){
	for(file in files){
		if(!str_detect(file, 'sig_genes') && !str_detect(file, 'png') && !str_detect(file, 'tiff')){
			pheno <- str_split(file, '/')[[1]][8]
			print(pheno)
			print(file)		
		if(str_detect(file, "ALL")){
        		model <- "ALL"
        		BP_Chrome <- load_bp_chrom('ALL', model_type)
      		} 
		if(str_detect(file, 'CAU')){
			model <- 'CAU'
			BP_Chrome <- load_bp_chrom('CAU', model_type)
		} 
		if(str_detect(file, 'CHN')){
			model <- 'CHN'
			BP_Chrome <- load_bp_chrom('CHN', model_type)
		} 
		if(str_detect(file, 'HIS')){
			model <- 'HIS'
			BP_Chrome <- load_bp_chrom('HIS', model_type)
		}
		if(str_detect(file, 'AFA')){
			model <- 'AFA'
			BP_Chrome <- load_bp_chrom('AFA', model_type)
		}

		S_Pred_file <- read.table(file, header = T,  sep = ',')      
      		S_Pred_file$GENE <- gsub("\\..*","",S_Pred_file$gene) 
         	S_Pred_file <- S_Pred_file[-c( 9)]
      		names(S_Pred_file)[names(S_Pred_file) == 'pvalue'] <- 'P'

      		GWAS <- merge(S_Pred_file, BP_Chrome, by = c('GENE', 'gene_name'))      
      		GWAS <- na.omit(GWAS)
      		names(GWAS)
		colnames(GWAS)[14] <- "CHR"
      	
      		GWAS <- GWAS %>%  #added by Jenny
        		transform(CHR = str_replace(CHR, "chr", "")) 
      		GWAS<- transform(GWAS, CHR=as.numeric(CHR)) 
		#elyse added this chunk to extract phenotype from my file name before it gets added in a column to the file
		if(str_detect(pheno, 'unfiltered')){
		pheno <- str_split(pheno, 'unfiltered_')[[1]][2]
		pheno <- str_replace(pheno, '.csv', '')
		} 
		if(str_detect(pheno, '0.05')){
		pheno <- str_split(pheno, '0.05_')[[1]][2]
		pheno <- str_replace(pheno, '.csv', '')
		} 
      		GWAS$Phenotype <- rep(pheno, nrow(GWAS))
      		GWAS$Model <- rep(model, nrow(GWAS))
      
      		threshold <- 5e-08
      		signif <- filter(GWAS,-log10(P) > -log10(5e-08))
      		print(signif)
      		signif <- unique(signif) #these significant hits aren't written out anywhere i don't think? separate bonferroni script for this 
		output <- rbind(output, GWAS)
#		if(dim(signif)[1] != 0){
#      			output <- rbind(output, GWAS)
#			return(GWAS)
#		}
		
		}
	}
	return(output)
}

make_manhattan <- function(output, model_type, num_pheno){
	axis.set <- output %>% 
  		group_by(CHR) %>% 
  		summarize(center = (max(BP) + min(BP)) / 2)
	ylim <- abs(floor(log10(min(output$P)))) + 2 
	if(ylim < 10){
		ylim <- 10
	}
	sig <- 5e-8

	# Prepare the dataset
	output <- output %>%  
 	 	# Compute chromosome size
 	 	group_by(CHR) %>% 
 	 	summarise(chr_len=max(BP)) %>% 
  
 	 	# Calculate cumulative position of each chromosome
 	 	mutate(tot=cumsum(chr_len)-chr_len) %>%
 	 	select(-chr_len) %>%
  
 	 	# Add this info to the initial dataset
  		left_join(output, ., by=c("CHR"="CHR")) %>%
  
  		# Add a cumulative position of each SNP
  		arrange(CHR, BP) %>%
  		mutate( BPcum=BP+tot) %>%
    
  		# Filter SNP to make the plot lighter
		filter(-log10(P)>0.5)

		# Prepare X axis
		axisdf <- output %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#		colorsss <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#276419", "#B15928", "#000000", "#8E0152")
		p <- ggplot(output, aes(x = BPcum, y = -log10(P), 
                          color = Phenotype, shape = Phenotype)) + facet_wrap(~Model, nrow = 5)+
  			geom_point(alpha = 0.75) +
  			geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  			scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) + 
  			scale_shape_manual(values=rep(c(15,16,17,18,4,5,6,7,8,9,10,15,16,17), 2)) + scale_color_manual(values = wes_palette("Zissou1", 28, type = "continuous")) +
  			scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  			scale_size_continuous(range = c(0.5,3)) +
  			labs(title = model_type, x = NULL, 
       			y = "-log10(p)") + 
  			theme_minimal() +
  			theme( legend.position = 'bottom', legend.title = element_text(size = 12), legend.text = element_text(size = 12), strip.text = element_text(size=13),  
    			panel.border = element_blank(),
    			panel.grid.major.x = element_blank(),
    			panel.grid.minor.x = element_blank(),
    			axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  			) + ggsave(paste(model_type, 'Manhattan_total_phenotype_5e8.tiff', sep = ''),device='tiff', width = 12, height = 8)

}


main <- function(){
	phenotypes <- c('LDL_cholesterol', 'Fasting_glucose', 'Fasting_insulin', 'Coffee', 'BMI', 'C-reactive', 'Chronic_kidney', 'Diabetes', 'Diastolic_blood', 'Hemoglobin', 'Glomerular','Mean_corpuscular_hemoglobin', 'HDL_cholesterol', 'Total_cholesterol', 'Triglyceride','Systolic_blood', 'Smoking', 'WBC', 'Height','Waist-hip50', 'Waist-hip51', 'Waist-hip52', 'End_renal', 'QRS_duration', 'QT_interval', 'PR_interval', 'Hypertension', 'Platelet')
	pops <- c('AFA', 'HIS', 'CAU', 'ALL', 'CHN')
	#model_types <- c('PCAIRbaseline_PAV_filtered_rho0.1_zpval0.05', 'PCAIRbaseline_PAV_filtered_unfiltered', 'PCAIRbaseline_rho0.1_zpval0.05', 'PCAIRbaseline_unfiltered')
	#must run one model at a time for reasons unknown. if all models given at once, the output ends up appending to the previous model's output, v bad, don't know why
	model_types <- c('PCAIRbaseline_rho0.1_zpval0.05')
	output <- data.frame()

	for(model_type in model_types){
		print(model_type)
		files <- get_files(model_type)
		print(files)
		output <- format_files(files, model_type, output)
		print(unique(output$Phenotype))
		write.csv(output, paste('/home/isabelle/topmed/proteome/PAGE/', model_type, '_sig_genes_all_pheno.csv', sep = ''), quote = F, row.names=F)
		#not rly significant genes, just all of them, idk man the name is what it is
		make_manhattan(output, model_type, len(unique(output$Phenotype)))		
	}
}

main()
