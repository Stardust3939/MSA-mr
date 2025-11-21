library(TwoSampleMR)
library(QTLMR)

library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)

library(mr.raps)
library(TwoSampleMR)
#library(plyr)
library(dplyr)
#library(fs) 
#library(ggplot2) 
#library(lubridate) 
library(ieugwasr)
library(plinkbinr) 
library(tidyr)
#library(progress)
library(data.table)
library(MRPRESSO)
library(parallel)
#library(foreach)
library(doParallel)
library(pbapply)
library(stringr)
library(cause)
library(vroom)
library(MungeSumstats)
#library(GenomicFiles)
#library(meta)
#library(readr)
#library(readxl)
#library(forestploter)
library(ldscr)
library(MRlap)

setwd("D:\\data\\MR_N_ukbppp")

pqtl_dat = data_cispQTL_ukb_ppp_exp(clump = 2)

pqtl_dat <- remove_MHC_data(dat = pqtl_dat,
                            chr_col = "chr.exposure",
                            pos_col = "pos.exposure",
                            MHC_start = 28477797,
                            MHC_end = 33448354)


for (file in files) {
  filename = basename(file)
 
  #dir.create(paste0("D:\\data\\MR_N_ukbppp\\", filename, "_MR"))
  
  #setwd(paste0("D:\\data\\MR_N_ukbppp\\", filename, "_MR"))
  #outcome_data <- read.table(paste0("D:\\data\\UKB_formated\\", filename, "_TwosampleMR.txt"), header = T, sep = "\t")
  #outcome_data <- read.csv(paste0("D:\\data\\UKB_formated\\", filename, "_TwosampleMR.txt"), header = T, sep = "\t")
  
  outcome_data <- modified_proxy_1000G(
    pqtl_dat$SNP,
    GWASfile = paste0("D:\\data\\UKB_formated\\", filename, "_TwosampleMR.txt"),
    proxies = T,
    rsq = 0.8,
    kb = 5000,
    nsnp = 5000,
    maf_threshold = 0.3,
    bfile_1000G = "D:\\data\\1kgv3\\EUR"
  )
  write.table(outcome_data, file = paste0(filename, "_outcome_data.txt"), sep = "\t", row.names = F)

}



for (file in files) {
  filename = basename(file)
  
  setwd("D:\\data\\MR_N_ukbppp")
  
  outcome_data <- read.table(paste0(filename, "_outcome_data.txt"), header = T, sep = "\t")
  
  dat <- harmonise_data(exposure_dat = pqtl_dat, outcome_dat = outcome_data)
  
  dir.create(paste0("D:\\data\\MR_N_ukbppp\\", filename, "_MR"))
  setwd(paste0("D:\\data\\MR_N_ukbppp\\", filename, "_MR"))
  
  res <- xQTL_mr(dat, FDR_method = "fdr", PVE = TRUE)
  
  write.table(res, file = paste0(filename, "_MR_results.txt"), sep = "\t", row.names = F)
}

files = list.dirs("D:\\data\\MR_N_ukbppp", full.names = T, recursive = F)
for (file in files) {
  setwd(file)
  filename = basename(file)
  res <- read.table(paste0(file,"/",filename, "_results.txt"), header = T, sep = "\t")
  xQTL_volcano_plot(
    res,
    breaks = c(0.1, 0.2, 0.3, 0.4, 0.5),
    scale_color = c("green", "grey", "red"),
    FDR_method = "fdr",
    save_plot = TRUE,
    pdf_name = "volcano_plot",
    width = 6,
    height = 5,
    save_path = "./"
  )
}

setwd("D:\\data\\MR_Net_ukbppp")

files <- list.files("D:\\data\\5775047", pattern = ".fastGWA", full.names = T)
for (file in files) {
  filename = basename(file)
  dir.create(paste0("D:\\data\\MR_Net_ukbppp\\", filename))
  setwd(paste0("D:\\data\\MR_Net_ukbppp\\", filename))
  
  outcome_dat <- read_outcome_data(
    filename = file,
    sep = "\t",
    #phenotype_col = "Phenotype",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "AF1",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "P",
    #units_col = "units",
    #ncase_col = "ncase",
    #ncontrol_col = "ncontrol",
    #samplesize_col = "samplesize",
    #gene_col = "gene",
    #id_col = "id",
    min_pval = 1e-200,
    log_pval = FALSE,
    chr_col = "CHR",
    pos_col = "POS"
  )
  
  write.table(outcome_dat, file = paste0(filename, "_outcome_data.txt"), sep = "\t", row.names = F)
  
  outcome_data <- modified_proxy_1000G(
    pqtl_dat$SNP,
    GWASfile = paste0("./", filename, "_outcome_data.txt"),
    proxies = T,
    rsq = 0.8,
    kb = 5000,
    nsnp = 5000,
    maf_threshold = 0.3,
    bfile_1000G = "D:\\data\\1kgv3\\EUR"
  )
  
  dat <- harmonise_data(exposure_dat = pqtl_dat, outcome_dat = outcome_data)
  res <- xQTL_mr(dat, FDR_method = "fdr", PVE = TRUE)
  write.table(res, file = paste0(filename, "_MR_results.txt"), sep = "\t", row.names = F)
  xQTL_volcano_plot(
    res,
    breaks = c(0.1, 0.2, 0.3, 0.4, 0.5),
    scale_color = c("green", "grey", "red"),
    FDR_method = "fdr",
    save_plot = TRUE,
    pdf_name = "volcano_plot",
    width = 6,
    height = 5,
    save_path = "./"
  )
}
