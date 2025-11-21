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


setwd("D:\\data\\MR_Vol_Meta")

# read exposure data:
msa <- read.table("D:\\MSA_TwosampleMR.txt", header = T, sep = "\t")

files <- list.files("D:\\data\\META_Vol", pattern = ".txt", full.names = T)

for (file in files) {
  tran_chrpos_from_SNP(
    GWASfile = file,
    snp_col = "MarkerName",
    chr_col = "chr",
    pos_col = "pos",
    build = 37,
    save_path = "./chrpos",
    save_name = paste0(basename(file), "_chrpos.txt")
  )
}

for (file in files) {
  # load outcome data:
  ukb <- read.table(file, header = T, sep = "\t")
  
  ukb <- read_outcome_data(
    file,
    snps = NULL,
    sep = "\t",
    #phenotype_col = "Phenotype",
    snp_col = "MarkerName",
    #beta_col = "BETA",
    #se_col = "SE",
    #eaf_col = "eaf",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    pval_col = "P-value",
    #units_col = "units",
    #ncase_col = "ncase",
    #ncontrol_col = "ncontrol",
    samplesize_col = "Weight",
    #gene_col = "gene",
    #id_col = "id",
    min_pval = 1e-200,
    log_pval = FALSE,
    chr_col = "CHR",
    pos_col = "BP"
  )
  
  # add se column to ukb:
  ukb$se <- sqrt(((ukb$beta.outcome)^2)/qchisq(ukb$pval.outcome,1,lower.tail=F))
  
  ukb <- get_eaf_from_1000G(ukb, "D:\\data", type = "outcome")
  
  # keep filename only:
  filename <- basename(file)
  
  # format outcome data:
  format_dat(
    ukb,
    type = "outcome",
    snps = NULL,
    header = TRUE,
    #phenotype_col = "Phenotype",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se",
    eaf_col = "eaf.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    #units_col = "units",
    #ncase_col = "n",
    #ncontrol_col = "ncontrol",
    #samplesize_col = "n",
    #gene_col = "gene",
    #id_col = "id",
    min_pval = 1e-200,
    #z_col = "z",
    #info_col = "info",
    chr_col = "chr.outcome",
    pos_col = "pos.outcome",
    log_pval = FALSE,
    Twosample_dat = TRUE,
    SMR_dat = FALSE,
    MTAG_dat = FALSE,
    METAL_dat = FALSE,
    GWAS_name = filename,
    save_path = "D:\\data\\UKB_formated"
  )
  #ukb_formated <- read.table(paste0("D:\\UKB_se_eaf\\", filename, "_TwosampleMR.txt"), header = T, sep = "\t")
  
}