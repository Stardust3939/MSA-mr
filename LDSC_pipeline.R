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

# get outcome file list:
files <- list.files("D:\\data\\UKB_treated", pattern = ".tsv", full.names = T)
# format all outcome files:
for (file in files) {
  # load outcome data:

  ukb <- read_outcome_data(
    file,
    snps = NULL,
    sep = "\t",
    #phenotype_col = "Phenotype",
    snp_col = "SNP",
    beta_col = "BETA",
    #se_col = "SE",
    #eaf_col = "eaf",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "P",
    #units_col = "units",
    #ncase_col = "ncase",
    #ncontrol_col = "ncontrol",
    samplesize_col = "NMISS",
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
    samplesize_col = "samplesize.outcome",
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
    save_path = "D:\\data\\UKB_formated2"
  )
  #ukb_formated <- read.table(paste0("D:\\UKB_se_eaf\\", filename, "_TwosampleMR.txt"), header = T, sep = "\t")

}

files = list.files("D:\\data\\UKB_formated2", pattern = ".txt", full.names = T)

for (file in files) {
  ldsc_rg <- LDSC_ldscr_rg(
    GWASfile = file,
    ancestry = "EUR",
    sample_prev = NA,
    population_prev = NA,
    ld,
    wld,
    n_blocks = 200,
    chisq_max = NA
  )

  ldsc_res <- LDSC_ldscr_h2(
    GWASfile = file,
    ancestry = "EUR",
    sample_prev = NA,
    population_prev = NA,
    ld,
    wld,
    n_blocks = 200,
    chisq_max = NA
  )
}

setwd("D:\\data\\LDSC")

filelist <- list.files("D:\\data\\UKB_formated2", pattern = ".txt", full.names = T)

for (file in filelist) {
  exposure <- read.table(file, header = T, sep = "\t")

  tryCatch({
    res <- LDSC_MRlap(
      exposure,
      exposure_name = "exposure",
      outcome,
      outcome_name = "outcome",
      ld = "D:\\data\\GenomicSEMFiles\\eur_w_ld_chr",
      hm3 = "D:\\data\\GenomicSEMFiles\\eur_w_ld_chr\\w_hm3.snplist",
      do_pruning = TRUE,
      user_SNPsToKeep = "",
      MR_threshold = 1e-5,
      MR_pruning_dist = 10000,
      MR_pruning_LD = 0.001,
      MR_reverse = 0.001,
      bfile_1000G = "D:\\data\\1kgv3\\EUR",
      save_logfiles = FALSE,
      verbose = TRUE
    )
  }, error = function(e) {
    print(paste0(file,"ERROR"))
  }, finally = {
    print(paste0(file,"OK"))
  })
}


exposure <- read.table(file, header = T, sep = "\t")
outcome <- read.table("D:\\data\\MSA_outcome.txt", header = T, sep = "\t")


res <- LDSC_MRlap(
  exposure,
  exposure_name = "exposure",
  outcome,
  outcome_name = "outcome",
  ld = "D:\\data\\GenomicSEMFiles\\eur_w_ld_chr",
  hm3 = "D:\\data\\GenomicSEMFiles\\eur_w_ld_chr\\w_hm3.snplist",
  do_pruning = TRUE,
  user_SNPsToKeep = "",
  MR_threshold = 1e-5,
  MR_pruning_dist = 10000,
  MR_pruning_LD = 0.001,
  MR_reverse = 0.001,
  bfile_1000G = "D:\\data\\1kgv3\\EUR",
  save_logfiles = FALSE,
  verbose = TRUE
)

res = MRlap(
  exposure = exposure,
  exposure_name = "exposure",
  outcome = outcome,
  outcome_name = "outcome",
  ld = "D:\\data\\GenomicSEMFiles\\eur_w_ld_chr",
  hm3 = "D:\\data\\GenomicSEMFiles\\eur_w_ld_chr\\w_hm3.snplist",
  MR_threshold = 1e-4,
  MR_pruning_LD = 0.05,
)

LDSC_py_sumstats(
  GWASfile = "D:\\data\\MSA_outcome.txt",
  GWAS_name = "MSA",
  N = NULL,
  merge_alleles = "D:\\data\\GenomicSEMFiles\\eur_w_ld_chr\\w_hm3.snplist",
  save_path = "./"
)

LDSC_py_sumstats(
  GWASfile = "D:\\data\\UKB_formated2\\ukb_roi_volume_may12_2019_phase1and2_pheno1_allchr_withA2.tsv_TwosampleMR.txt",
  GWAS_name = "outcome",
  N = NULL,
  merge_alleles = "D:\\data\\GenomicSEMFiles\\eur_w_ld_chr\\w_hm3.snplist",
  save_path = "./"
)

LDSC_py_rg(
  test_help = FALSE,
  Sumstatsfile = c("./outcome.sumstats.gz","./MSA.sumstats.gz"),
  ref_ld_chr = "D:\\data\\GenomicSEMFiles\\eur_w_ld_chr",
  w_ld_chr = "D:\\data\\GenomicSEMFiles\\eur_w_ld_chr",
  opt_arguments = NULL,
  save_name = "res",
  save_path = "./"
)

files = list.files("D:\\data\\UKB_formated2", pattern = ".txt", full.names = T)

for (file in files) {
  filename = basename(file)
  LDSC_py_sumstats(
    GWASfile = file,
    GWAS_name = filename,
    N = NULL,
    merge_alleles = "D:\\data\\GenomicSEMFiles\\eur_w_ld_chr\\w_hm3.snplist",
    save_path = "./"
  )
}
files = list.files("D:\\data\\LDSC", pattern = ".gz", full.names = T)
for (file in files) {
  filename = basename(file)
  LDSC_py_rg(
    test_help = FALSE,
    Sumstatsfile = c("D:\\data\\MSA.sumstats.gz",file),
    ref_ld_chr = "D:\\data\\GenomicSEMFiles\\eur_w_ld_chr",
    w_ld_chr = "D:\\data\\GenomicSEMFiles\\eur_w_ld_chr",
    opt_arguments = NULL,
    save_name = filename,
    save_path = "./"
  )
}
setwd("D:\\data\\Net_LDSC")
files <- list.files("D:\\data\\5775047", pattern = ".fastGWA", full.names = T)
for (file in files) {
  filename = basename(file)

  exposure_dat <- read_exposure_data(
    filename = file,
    clump = FALSE,
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
    samplesize_col = "N",
    #gene_col = "gene",
    #id_col = "id",
    min_pval = 1e-200,
    log_pval = FALSE,
    chr_col = "CHR",
    pos_col = "POS"
  )
  write.table(exposure_dat, paste0("D:\\data\\net_formated\\", filename, "_TwosampleMR.txt"), sep = "\t", row.names = F)
}
files = list.files("D:\\data\\net_formated", pattern = ".txt", full.names = T)
for (file in files) {
  filename = basename(file)
  LDSC_py_sumstats(
    GWASfile = file,
    GWAS_name = filename,
    N = NULL,
    merge_alleles = "D:\\data\\GenomicSEMFiles\\eur_w_ld_chr\\w_hm3.snplist",
    save_path = "./"
  )
}
files = list.files("D:\\data\\Net_LDSC", pattern = ".gz", full.names = T)
for (file in files) {
  filename = basename(file)
  LDSC_py_rg(
    test_help = FALSE,
    Sumstatsfile = c("D:\\data\\MSA.sumstats.gz",file),
    ref_ld_chr = "D:\\data\\GenomicSEMFiles\\eur_w_ld_chr",
    w_ld_chr = "D:\\data\\GenomicSEMFiles\\eur_w_ld_chr",
    opt_arguments = NULL,
    save_name = filename,
    save_path = "./"
  )
}
