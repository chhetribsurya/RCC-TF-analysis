



################
##Extract summary statistics for each RCC SNP-protein combination
################

beta <- sd <- af <- NULL

for(prot in prots){
  
  
  tx <- read.table(paste0("/data/Dutta_lab/decode/ukb_pval/renal_cell_carcinoma/",prot,".txt"),header = T)
  beta <- cbind(beta,tx$beta)
  sd <- cbind(sd,tx$sd)
  af <- cbind(af,tx$af)
  
  print(prot)
  
}


#################


tfs = read.table("/data/Dutta_lab/decode/rcc_TF_paper/TFBS-SNP-Matrix/rcc_tfbs_matrix_10000/Renal_cell_carcinoma_output_FULLmatrix.txt",header = T)
tfs = tfs[match(t.loci$rsid, tfs$snpID),]

indx <- which(colSums(tfs[,-1]) == 0)
tfs <- tfs[,-(indx+1)]

list.indx <- 1
fin.pv <- list()

for(j in 1:dim(tfs)[2]){
  
  wts = wts.orig*tfs[,j]
  
  pv.skat <- pv.burden <- NULL
  
  for(i in 1:length(t2)){
    
    bx = beta[,i]; sdx = sd[,i]; afx = af[,i];
    
    
    qx =  kernel_stats(bx,sdx,nsamp = 34557,kernel = diag(length(bx)),weights = wts)
    A1 = get.w.mat(afx,n.samp = 34557,is.burden = FALSE,weights = wts)
    pv.skat <- c(pv.skat, SKAT:::Get_PValue.Lambda(lambda = eigen(A1)$values,Q = qx)$p.value)
    
    qx =  kernel_stats(bx,sdx,nsamp = 34557,kernel = matrix(1,ncol = length(bx),nrow = length(bx)),weights = wts)
    A1 = get.w.mat(afx,n.samp = 34557,is.burden = TRUE,weights = wts)
    pv.burden <- c(pv.burden, SKAT:::Get_PValue.Lambda(lambda = eigen(A1)$values,Q = qx)$p.value)
    
    print(i)
    
  }
  names(pv.skat) = prots
  names(pv.burden) = prots
  
  fin.pv[[list.indx]]= data.frame(pv.burden,pv.skat)
  print(paste0("Done TF ",colnames(tfs)[j]))
  list.indx <- list.indx + 1
}

names(fin.pv) <- colnames(tfs)[-1]


trans.list <- list()

for(i in 1:dim(t.loci)[1]){  
  
  ss.sub = prot.info %>% filter(chr == t.loci$chromosome[i])
  endx <- abs(ss.sub$gene_end - t.loci$base_pair_location[i])
  stx <- abs(ss.sub$gene_start - t.loci$base_pair_location[i])
  px <- union(which(stx < 1000000), which(endx < 1000000))  
  trans.list[[i]] <- ss.sub$OlinkID[px]
  print(i)  
  
}


#293_output_FULLmatrix.txt      786-O_output_FULLmatrix.txt  Kidney_Cortex_output_FULLmatrix.txt  Renal_cell_carcinoma_output_FULLmatrix.txt
#786-M1A_output_FULLmatrix.txt  G-401_output_FULLmatrix.txt  RCC4_output_FULLmatrix.txt           TTC1240_output_FULLmatrix.txt


##################



t.loci <- read.table("/data/Dutta_lab/decode/RCC_loci.txt",header = T)


# indx <- which(colSums(tfs[,-1]) == 0)
# tfs <- tfs[,!indx]


setwd("/data/Dutta_lab/decode/ukb_ppp/EUR/")
prot.all <- system("ls -d *",intern = T)
prot = gsub(pattern="/",replacement="",prot.all)
t2 <- prot

# betas <- read.table("/data/Dutta_lab/wangky/UKB_project/All_ptns_RCC_loci/BETA",header = T)
# sd <- read.table("/data/Dutta_lab/wangky/UKB_project/All_ptns_RCC_loci/SE",header = T)
# af <- read.table("/data/Dutta_lab/wangky/UKB_project/All_ptns_RCC_loci/MAF",header = T)


beta = sd = af = NULL
for(i in 1:length(t2)){
  
  infile <- paste0("/data/Dutta_lab/decode/ukb_pval/renal_cell_carcinoma/",t2[i],".txt");
  tx <- read.table(infile,header = T);
  beta <- cbind(beta,tx[,1]); sd <- cbind(sd,tx[,2]); af <- cbind(af,tx[,3])
  print(i)
}

colnames(beta) <- colnames(sd) <- colnames(af) <- t2
wts <- t.loci$beta

wts.orig <- wts



library(tidyverse)
library(SKAT)
source("/data/Dutta_lab/decode/codes/group_test_skat.R")


infile <- "/data/Dutta_lab/decode/rcc_TF_paper/TFBS-SNP-Matrix/rcc_tfbs_matrix_10000/293_output_FULLmatrix.txt"
outfile <- "/data/Dutta_lab/decode/rcc_TF_paper/results_TF_subPRS/293_output_10000.rds"

tfs = read.table(infile,header = T)
tfs = tfs[match(t.loci$rsid, tfs$snpID),]

list.indx <- 1
fin.pv <- list()


for(j in 1:dim(tfs)[2]){
  
  wts = wts.orig*tfs[,j]
  
  pv.skat <- pv.burden <- NULL
  
  for(i in 1:length(t2)){
    
    bx = beta[,i]; sdx = sd[,i]; afx = af[,i];
    
    
    qx =  kernel_stats(bx,sdx,nsamp = 34557,kernel = diag(length(bx)),weights = wts)
    A1 = get.w.mat(afx,n.samp = 34557,is.burden = FALSE,weights = wts)
    pv.skat <- c(pv.skat, SKAT:::Get_PValue.Lambda(lambda = eigen(A1)$values,Q = qx)$p.value)
    
    qx =  kernel_stats(bx,sdx,nsamp = 34557,kernel = matrix(1,ncol = length(bx),nrow = length(bx)),weights = wts)
    A1 = get.w.mat(afx,n.samp = 34557,is.burden = TRUE,weights = wts)
    pv.burden <- c(pv.burden, SKAT:::Get_PValue.Lambda(lambda = eigen(A1)$values,Q = qx)$p.value)
    
    print(i)
    
  }
  names(pv.skat) = prot
  names(pv.burden) = prot
  
  fin.pv[[list.indx]]= data.frame(pv.burden,pv.skat)
  print(paste0("Done TF ",colnames(tfs)[j]))
  list.indx <- list.indx + 1
}

names(fin.pv) <- colnames(tfs)[-1]

saveRDS(fin.pv,outfile)

