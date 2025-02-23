#!/usr/bin/env R

library(data.table)
library(dplyr)
library(tidyverse)
library(lme4)

#pass arguments:
sysargv <- commandArgs(TRUE)
bedfile_1 <- sysargv[1]
bedfile_2 <- sysargv[2]
workdir <- sysargv[3]
outputdir <- sysargv[4]

#assign argvariables
tf_bedfile_1 <- bedfile_1
tf_bedfile_2 <- bedfile_2
work_dir <- workdir
output_dir <- outputdir

# label file and dirs
rcc_gwas <- c("RCC", "CC2", "PRCC")
#rcc_gwas <- c("RCC")
#rcc_gwas <- c("RCC", "CC2")

#work_dir <- "/Users/suryachhetri/Dropbox/for_diptavo/RCC_TF_Project"
#work_dir <- "/scratch16/abattle4/surya/datasets/for_diptavo/RCC_TF_Project"
gwas_dir <- file.path(work_dir, "GWAS_RCC")


# create empty list
zscore_based <- list()
chisqr_based <- list()
zscore_beta_vals <- list()
zscore_se_vals <- list()
chisqr_beta_vals <- list()
chisqr_se_vals <- list()
#gwas_rcc <- list()
#ldblock <- list()

tf_name_1 <- strsplit(basename(tf_bedfile_1), "\\.")[[1]][1]
tf_name_2 <- strsplit(basename(tf_bedfile_2), "\\.")[[1]][1]
cat("\n\n")
print(paste("Processing TF : ", tf_name_1))
print(paste("Processing TF : ", tf_name_2))

for (gwas_type in rcc_gwas) {

    snpbedformat_file <- file.path(gwas_dir, paste0("MultiAnc_", gwas_type, ".BedToolformat.bed"))
    snpbedsorted_file <- file.path(gwas_dir, paste0("MultiAnc_", gwas_type, ".BedToolformat.srt.bed"))

    # create dirs
    #out_dir <- file.path(work_dir, "HEK293_BATCH_output")
    #out_dir <- file.path(work_dir, "HEK293_BATCH_output", "tf_cooccurence")
    #out_dir <- file.path(work_dir, "test_HEK293_BATCH_output", "tf_cooccurence")
    out_dir <- file.path(output_dir, "HEK293_BATCH_output_expanded", "tf_cooccurence")
    sort_dir <- file.path(out_dir, "sort_intersect", paste0(tf_name_1,"-",tf_name_2))
    if (!dir.exists(out_dir)){dir.create(out_dir, recursive = TRUE)}
    if (!dir.exists(sort_dir)){dir.create(sort_dir, recursive = TRUE)}

    # read gwas files and generate snp bed
    gwas_df <- fread(file.path(work_dir, "GWAS_RCC", paste0("MultiAnc_", gwas_type, ".tsv")))
    snp_bed <- data.frame(paste0("chr", gwas_df$CHR), gwas_df$BP, as.integer(gwas_df$BP+1), gwas_df$SNP)

    if (!file.exists(snpbedformat_file)){
        write.table(snp_bed, snpbedformat_file, col.names=F, row.names=F, quote = F, sep = "\t")
    }

    # bash command and bedtools sort for snp overlap with TFs
    # activate conda environment "bedtools" and run command
    CMD_sort <- paste("bedtools", "sort", "-i", snpbedformat_file, ">", snpbedsorted_file)
    if (!file.exists(snpbedsorted_file)){
        system(paste("conda run -n bedtools", CMD_sort))
    }
    #tfbed_dir <- file.path(work_dir, "HEK293_bed")
    #tf_bed <- "ATF2.bed"
    #bedlist <- list.files(tfbed_dir, pattern=".bed")

    cat("\n")
    print(paste("#Processing GWAS TYPE: ", gwas_type))
    cat("\n")
    
    #########################
    #TF cooccurence analysis
    #TF1 status compute
    #########################

    #tf_bedfile <- file.path(tfbed_dir, tf_bed)
    tf_bedsort_file_1 <- file.path(sort_dir, paste0(tf_name_1, ".srt.bed"))
    snp_tf_file_1 <- file.path(sort_dir, paste0(tf_name_1, "_SNP_intsct.bed"))
    
    # bash command lines
    CMD_bedsort <- paste("bedtools", "sort", "-i", tf_bedfile_1, ">", tf_bedsort_file_1)
    CMD_intsect <- paste("bedtools", "intersect", "-a", snpbedsorted_file, "-b", tf_bedsort_file_1, ">", snp_tf_file_1)

    # activate conda environment "bedtools" and run command
    #if (!file.exists(tf_bedsort_file)){
    #    system(paste("conda run -n bedtools", CMD_bedsort))
    #}
    
    system(paste("conda run -n bedtools", CMD_bedsort))
    
    #if (!file.exists(snp_tf_file)){
    #    system(paste("conda run -n bedtools", CMD_intsect))
    #}

    system(paste("conda run -n bedtools", CMD_intsect))
    
    tf_bed_1 <- read.table(snp_tf_file_1)
    tf_status_1 <- ifelse(gwas_df$SNP %in% tf_bed_1$V4, yes = 1, no = 0)

    #########################
    #TF cooccurence analysis
    #TF2 status compute
    #########################

    #tf_bedfile <- file.path(tfbed_dir, tf_bed)
    tf_bedsort_file_2 <- file.path(sort_dir, paste0(tf_name_2, ".srt.bed"))
    snp_tf_file_2 <- file.path(sort_dir, paste0(tf_name_2, "_SNP_intsct.bed"))
    
    # bash command lines
    CMD_bedsort <- paste("bedtools", "sort", "-i", tf_bedfile_2, ">", tf_bedsort_file_2)
    CMD_intsect <- paste("bedtools", "intersect", "-a", snpbedsorted_file, "-b", tf_bedsort_file_2, ">", snp_tf_file_2)

    # activate conda environment "bedtools" and run command
    #if (!file.exists(tf_bedsort_file)){
    #    system(paste("conda run -n bedtools", CMD_bedsort))
    #}
    
    system(paste("conda run -n bedtools", CMD_bedsort))
    
    #if (!file.exists(snp_tf_file)){
    #    system(paste("conda run -n bedtools", CMD_intsect))
    #}

    system(paste("conda run -n bedtools", CMD_intsect))
    
    tf_bed_2 <- read.table(snp_tf_file_2)
    tf_status_2 <- ifelse(gwas_df$SNP %in% tf_bed_2$V4, yes = 1, no = 0)
    # accounting LD block effects on TFs association
    #kb <- 1000000 ## also we can try 500000, 1000000, 250000

    #ldblock_list <- c(2500000, 2000000, 1000000, 500000, 250000)
    ldblock_list <- c(2500000, 2000000, 1000000)
    ldblock_list <- c(2500000)

    for (kb in ldblock_list){
        
        print(paste("Processing LDBlock: ", toString(kb)))
        gwas_df$grp <- paste0(gwas_df$CHR, '_', floor(gwas_df$BP/kb))

        # our approach
        # zscore <- qnorm(gwas_df$p, lower.tail=F)
        pval <- gwas_df$p
        pval[which(pval == 1)] <- 0.99 #set it to 0.99 to counter infinity error on y
        zscore <- qnorm(pval, lower.tail=F)
        lm_res_1 <- lmer(I(zscore) ~ tf_status_1 + tf_status_2 + tf_status_1*tf_status_2 + (1|gwas_df$grp), control = lmerControl(calc.derivs = FALSE))
        lm_res_2 <- lmer(I(zscore) ~ tf_status_1*tf_status_2 + (1|gwas_df$grp), control = lmerControl(calc.derivs = FALSE))
        pvalue_1 <- pnorm(abs(summary(lm_res_1)$coef[4,3]), lower.tail=F)
        pvalue_2 <- pnorm(abs(summary(lm_res_2)$coef[2,3]), lower.tail=F)

        betas_1 <- summary(lm_res_1)$coefficients[4,2]
        betas_2 <- summary(lm_res_2)$coefficients[2,2]
        
        se_1 <- summary(lm_res_1)$coefficients[4,3]
        se_2 <- summary(lm_res_2)$coefficients[2,3]

        # append pvalues to pvalue list
        annot_1 <- paste(tf_name_1, tf_name_2, gwas_type, toString(kb), "TFindepEffectsAdj", sep="-")
        annot_2 <- paste(tf_name_1, tf_name_2, gwas_type, toString(kb), "TFinteractionOnly", sep="-")
        #annot <- paste(tf_name, gwas_type, toString(kb), sep="-")
        zscore_based[[annot_1]] <- pvalue_1
        zscore_based[[annot_2]] <- pvalue_2
        
        zscore_beta_vals[[annot_1]] <- betas_1
        zscore_beta_vals[[annot_2]] <- betas_2
        
        zscore_se_vals[[annot_1]] <- se_1
        zscore_se_vals[[annot_2]] <- se_2
        #gwas_rcc[[gwas_type]] <- gwas_type
        #ldblock[[toString(kb)]] <- toString(kb)

        # BRCA tf overlap paper approach
        #chisq <- (gwas_df$beta/gwas_df$standard_error)^2
        chisq <-  (gwas_df$Z)^2
        lm_res1 <- lmer(I(chisq) ~ tf_status_1 + tf_status_2 + tf_status_1*tf_status_2 + (1|gwas_df$grp), control = lmerControl(calc.derivs = FALSE))
        lm_res2 <- lmer(I(chisq) ~ tf_status_1*tf_status_2 + (1|gwas_df$grp), control = lmerControl(calc.derivs = FALSE))
        pvalue1 <- pnorm(abs(summary(lm_res1)$coef[4,3]), lower.tail=F)
        pvalue2 <- pnorm(abs(summary(lm_res2)$coef[2,3]), lower.tail=F)
            
        #betas_1 <- summary(lm_res_1)$coefficients[4,2]
        #betas_2 <- summary(lm_res_2)$coefficients[2,2]
        
        #se_1 <- summary(lm_res_1)$coefficients[4,3]
        #se_2 <- summary(lm_res_2)$coefficients[2,3]

        # append pvalues to pvalue list
        annot1 <- paste(tf_name_1, tf_name_2, gwas_type, toString(kb), "TFindepEffectsAdj", sep="-")
        annot2 <- paste(tf_name_1, tf_name_2, gwas_type, toString(kb), "TFinteractionOnly", sep="-")

        chisqr_based[[annot1]] <- pvalue1
        chisqr_based[[annot2]] <- pvalue2
            
        #chisqr_beta_vals[[annot_1]] <- betas_1
        #chisqr_beta_vals[[annot_2]] <- betas_2
        
        #chisqr_se_vals[[annot_1]] <- se_1
        #chisqr_se_vals[[annot_2]] <- se_2
        
        #task status update
        print(paste("#Task Completed Annotation: ", annot1))
        print(paste("#Task Completed Annotation: ", annot2))


    }

    print(paste("#Task Completed GWAS TYPE: ", gwas_type))
}


# combine list to dataframe
zscore_method <- unlist(zscore_based)
chisqr_method <- unlist(chisqr_based)

zscore_betas <- unlist(zscore_beta_vals)
#chisqr_betas <- unlist(chisqr_beta_vals)

zscore_se <- unlist(zscore_se_vals)
#chisqr_se <- unlist(chisqr_se_vals)

#gwas_rcc_type <- unlist(gwas_rcc)
#ldblock_type <- unlist(ldblock)

#cbind_df <- data.frame(zscore_method, chisqr_method)
#cbind_df <- data.frame(zscore_method, chisqr_method, zscore_betas, chisqr_betas, zscore_se, chisqr_se)
cbind_df <- data.frame(zscore_method, chisqr_method, zscore_betas, zscore_se)
#cbind_df <- data.frame(zscore_method, chisqr_method, gwas_rcc_type, ldblock_type)

final_df <- rownames_to_column(cbind_df, "TF_name")
#final_df_0.01 <- final_df %>% filter(zscore_method < 0.01)
#final_df <- cbind(zscore_list, chisqr_list)

# write SNP TF assoc file
annot <- paste(tf_name_1, tf_name_2, sep="-")
fwrite(final_df, file.path(out_dir, paste0(annot, ".txt")), sep="\t", col.names=F, row.names=F, quote=F)
#fwrite(final_df_0.01, file.path(out_dir, "SNP_TF_Assoc_final.0.01.txt"), sep="\t", col.names=T, row.names=T)

print("Final Task Completed...")

