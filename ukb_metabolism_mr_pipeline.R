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
#library(cause)
library(vroom)
library(MungeSumstats)
#library(GenomicFiles)
#library(meta)
#library(readr)
#library(readxl)
#library(forestploter)
library(ldscr)
library(MRlap)
# set work dir I:\\data
setwd("I:\\data\\ukb")

folders <- list.dirs("I:\\data\\ukb_metabolism", recursive = FALSE)

for (folder in folders) {
  setwd(folder)
  foldername = basename(folder)
  folder <- file.path(folder, "harmonised")
  files = list.files(folder, pattern = "*.tsv.gz", full.names = TRUE)
  files <- list.files(folder, pattern = "*.tsv.gz", full.names = TRUE)
  
  dat_Catalog <- format_data_Catalog(GWASfile = files[1],
                                     GWAS_name = basename(files[1]),
                                     samplesize = 450015, 
                                     get_SNP = FALSE,
                                     build = 38,
                                     save_path = "./outdata",
                                     type = "outcome",
                                     snp_col = "rsid",
                                     beta_col = "beta",
                                     se_col = "standard_error",
                                     eaf_col = "effect_allele_frequency",
                                     effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",
                                     pval_col = "p_value",
                                     chr_col = "chromosome",
                                     pos_col = "base_pair_location",
                                     Twosample_dat = TRUE,
                                     SMR_dat = FALSE,
                                     MTAG_dat = FALSE,
                                     METAL_dat = FALSE)
}

msaexposure <- read_exposure_data(
  "D:\\MSA_TwosampleMR.txt",
  sep = "\t",
  #phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  samplesize_col = "samplesize.exposure",
  #gene_col = "gene",
  #id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE,
  chr_col = "chr.exposure",
  pos_col = "pos.exposure"
)
exposure_dat <- msaexposure[msaexposure$pval.exposure <= 1e-4,]
colnames(exposure_dat)[colnames(exposure_dat) == "SNP"] <- "rsid"
colnames(exposure_dat)[colnames(exposure_dat) == "pval.exposure"] <- "pval"
colnames(exposure_dat)[colnames(exposure_dat) == "id.exposure"] <- "trait_id"


# msa_clump <- clump_data(msa, clump_r2 = 0.01, clump_kb = 5000, pop = "EUR")



exposure_clump <- ieugwasr::ld_clump(dat = exposure_dat, clump_r2 = 0.01, clump_kb = 5000, clump_p = 1, plink_bin = "C:/Users/Jerry/AppData/Local/R/win-library/4.4/plinkbinr/bin/plink_Windows.exe", bfile = "D:\\data\\1kgv3\\EUR")

colnames(exposure_clump)[colnames(exposure_clump) == "rsid"] <- "SNP"
colnames(exposure_clump)[colnames(exposure_clump) == "pval"] <- "pval.exposure"
colnames(exposure_clump)[colnames(exposure_clump) == "trait_id"] <- "id.exposure"





for (folder in folders[1:2]) {
  setwd(folder)
  foldername = basename(folder)
  folder <- file.path(folder, "outdata")
  files = list.files(folder, pattern = "*.txt", full.names = TRUE)
  outcome <- read_outcome_data(
    files,
    snps = NULL,
    sep = "\t",
    #phenotype_col = "Phenotype",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    eaf_col = "eaf.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    #units_col = "units",
    #ncase_col = "ncase",
    #ncontrol_col = "ncontrol",
    samplesize_col = "samplesize.outcome",
    #gene_col = "gene",
    #id_col = "id",
    min_pval = 1e-200,
    log_pval = FALSE,
    chr_col = "chr.outcome",
    pos_col = "pos.outcome"
  )
  filename = basename(files)
  dat <- harmonise_data(exposure_dat = exposure_clump, outcome_dat = outcome)
  res <- mr(dat)
  #export res to tsv in working directory:
  wd = getwd()
  write.table(res, file = paste0(wd, "\\", filename, "_MR_result.tsv"), sep = "\t", row.names = F)
  
  het <- mr_heterogeneity(dat)
  #export het to tsv in working directory:
  write.table(het, file = paste0(wd, "\\", filename, "_MR_heterogeneity.tsv"), sep = "\t", row.names = F)
  
  ple <- mr_pleiotropy_test(dat)
  #export ple to tsv in working directory:
  write.table(ple, file = paste0(wd, "\\", filename, "_MR_pleiotropy.tsv"), sep = "\t", row.names = F)
  
  #p1 = mr_scatter_plot(res,dat)
  png(file = "scatter_plot.png", width = 800, height = 800)
  p1 = mr_scatter_plot(res,dat)
  print(p1)
  dev.off()
  
  result_single <- mr_singlesnp(dat)
  #p2 = mr_forest_plot(result_single)
  png(file = "forest_plot.png", width = 800, height = 800)
  p2 = mr_forest_plot(result_single)
  print(p2)
  dev.off()
  
  result_loo <- mr_leaveoneout(dat)
  #p3 = mr_leaveoneout_plot(result_loo)
  png(file = "leaveoneout_plot.png", width = 800, height = 800)
  p3 = mr_leaveoneout_plot(result_loo)
  print(p3)
  dev.off()
  
  #p4=mr_funnel_plot(result_single)
  png(file = "funnel_plot.png", width = 800, height = 800)
  p4 = mr_funnel_plot(result_single)
  print(p4)
  dev.off()
  
  Visualizing_MR_forest(res,save_plot = TRUE,plot_pdf = "forest_plot")
}


msaoutcome <- read_outcome_data(
  "D:\\MSA_TwosampleMR.txt",
  snps = NULL,
  sep = "\t",
  #phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  samplesize_col = "samplesize.exposure",
  #gene_col = "gene",
  #id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE,
  chr_col = "chr.exposure",
  pos_col = "pos.exposure"
)


for (folder in folders) {
  setwd(folder)
  foldername = basename(folder)
  folder <- file.path(folder, "outdata")
  files = list.files(folder, pattern = "*.txt", full.names = TRUE)
  exposure_dat <- read_exposure_data(
    files,
    sep = "\t",
    #phenotype_col = "Phenotype",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    eaf_col = "eaf.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    #units_col = "units",
    #ncase_col = "ncase",
    #ncontrol_col = "ncontrol",
    samplesize_col = "samplesize.outcome",
    #gene_col = "gene",
    #id_col = "id",
    min_pval = 1e-200,
    log_pval = FALSE,
    chr_col = "chr.outcome",
    pos_col = "pos.outcome"
  )
  exposure_dat <- exposure_dat[exposure_dat$pval.exposure <= 1e-4,]
  colnames(exposure_dat)[colnames(exposure_dat) == "SNP"] <- "rsid"
  colnames(exposure_dat)[colnames(exposure_dat) == "pval.exposure"] <- "pval"
  colnames(exposure_dat)[colnames(exposure_dat) == "id.exposure"] <- "trait_id"
  
  
  # msa_clump <- clump_data(msa, clump_r2 = 0.01, clump_kb = 5000, pop = "EUR")
  
  
  
  exposure_clump <- ieugwasr::ld_clump(dat = exposure_dat, clump_r2 = 0.01, clump_kb = 5000, clump_p = 1, plink_bin = "C:/Users/Jerry/AppData/Local/R/win-library/4.4/plinkbinr/bin/plink_Windows.exe", bfile = "D:\\data\\1kgv3\\EUR")
  
  colnames(exposure_clump)[colnames(exposure_clump) == "rsid"] <- "SNP"
  colnames(exposure_clump)[colnames(exposure_clump) == "pval"] <- "pval.exposure"
  colnames(exposure_clump)[colnames(exposure_clump) == "trait_id"] <- "id.exposure"
  filename = basename(files)
  dat <- harmonise_data(exposure_dat = exposure_clump, outcome_dat = msaoutcome)
  res <- mr(dat)
  #export res to tsv in working directory:
  wd = getwd()
  dir.create(paste0(wd, "\\", "MSAoutcome"))
  setwd(paste0(wd, "\\", "MSAoutcome"))
  wd = getwd()
  write.table(res, file = paste0(wd, "\\", filename, "_MR_result.tsv"), sep = "\t", row.names = F)
  
  het <- mr_heterogeneity(dat)
  #export het to tsv in working directory:
  write.table(het, file = paste0(wd, "\\", filename, "_MR_heterogeneity.tsv"), sep = "\t", row.names = F)
  
  ple <- mr_pleiotropy_test(dat)
  #export ple to tsv in working directory:
  write.table(ple, file = paste0(wd, "\\", filename, "_MR_pleiotropy.tsv"), sep = "\t", row.names = F)
  
  #p1 = mr_scatter_plot(res,dat)
  png(file = "scatter_plot.png", width = 800, height = 800)
  p1 = mr_scatter_plot(res,dat)
  print(p1)
  dev.off()
  
  result_single <- mr_singlesnp(dat)
  #p2 = mr_forest_plot(result_single)
  png(file = "forest_plot.png", width = 800, height = 800)
  p2 = mr_forest_plot(result_single)
  print(p2)
  dev.off()
  
  result_loo <- mr_leaveoneout(dat)
  #p3 = mr_leaveoneout_plot(result_loo)
  png(file = "leaveoneout_plot.png", width = 800, height = 800)
  p3 = mr_leaveoneout_plot(result_loo)
  print(p3)
  dev.off()
  
  #p4=mr_funnel_plot(result_single)
  png(file = "funnel_plot.png", width = 800, height = 800)
  p4 = mr_funnel_plot(result_single)
  print(p4)
  dev.off()
  
  Visualizing_MR_forest(res,save_plot = TRUE,plot_pdf = "forest_plot")
}
