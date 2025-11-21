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

setwd("D:\\data\\MR_Vol")

# load exposure data:
msa <- read.table("D:\\data\\GCST90406925\\GCST90406925_hg19.h.txt", header = T, sep = "\t")

# format exposure data:

format_dat(
  msa,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  #units_col = "units",
  #ncase_col = "n",
  #ncontrol_col = "ncontrol",
  samplesize_col = "n",
  #gene_col = "gene",
  #id_col = "id",
  min_pval = 1e-200,
  #z_col = "z",
  #info_col = "info",
  chr_col = "chromosome",
  pos_col = "base_pair_location",
  log_pval = FALSE,
  Twosample_dat = TRUE,
  SMR_dat = FALSE,
  MTAG_dat = FALSE,
  METAL_dat = FALSE,
  GWAS_name = "MSA",
  save_path = "D:\\"
)

# read exposure data:
msa <- read.table("D:\\MSA_TwosampleMR.txt", header = T, sep = "\t")

# get outcome file list:
files <- list.files("D:\\data\\UKB_treated", pattern = ".tsv", full.names = T)
# format all outcome files:
for (file in files) {
  # load outcome data:
  ukb <- read.table(file, header = T, sep = "\t")
  
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

setwd("D:\\data\\MR_Vol")

for (file in files) {
  filename = basename(file)
  #create new folder
  
  dir.create(paste0("D:\\data\\MR_Vol\\", filename, "_MR"))
}

#rename columns SNP to rsid in dataframe msa:
colnames(msa)[colnames(msa) == "SNP"] <- "rsid"
colnames(msa)[colnames(msa) == "pval.exposure"] <- "pval"
colnames(msa)[colnames(msa) == "id.exposure"] <- "trait_id"


msa_clump <- clump_data(msa, clump_r2 = 0.01, clump_kb = 5000, pop = "EUR")



msa_clump <- ieugwasr::ld_clump(dat = msa, clump_r2 = 0.01, clump_kb = 5000, clump_p = 1, plink_bin = "C:/Users/Jerry/AppData/Local/R/win-library/4.4/plinkbinr/bin/plink_Windows.exe", bfile = "D:\\data\\1kgv3\\EUR")

colnames(msa_clump)[colnames(msa_clump) == "rsid"] <- "SNP"
colnames(msa_clump)[colnames(msa_clump) == "pval"] <- "pval.exposure"
colnames(msa_clump)[colnames(msa_clump) == "trait_id"] <- "id.exposure"

# get pval <= 0.001 in msa_clump:
msa_clump_pval <- msa_clump[msa_clump$pval.exposure <= 1e-4,]


for (file in files) {
  filename = basename(file)
  setwd(paste0("D:\\data\\MR_Vol\\", filename, "_MR"))
  #outcome_data <- read.table(paste0("D:\\data\\UKB_formated\\", filename, "_TwosampleMR.txt"), header = T, sep = "\t")
  outcome_data <- read.csv(paste0("D:\\data\\UKB_formated\\", filename, "_TwosampleMR.txt"), header = T, sep = "\t")
  
  dat <- harmonise_data(exposure_dat = msa_clump_pval, outcome_dat = outcome_data)
  
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

setwd(paste0("D:\\data\\MR_Vol\\", filename, "_MR"))

outcome_data <- read.table(paste0("D:\\data\\UKB_formated\\", filename, "_TwosampleMR.txt"), header = T, sep = "\t")

dat <- harmonise_data(exposure_dat = msa_clump_pval, outcome_dat = outcome_data)

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


for (file in files) {
  filename = basename(file)
  setwd(paste0("D:\\data\\MR_Vol\\", filename, "_MR"))
  res = read.table(paste0("D:\\data\\MR_Vol\\", filename, "_MR\\", filename, "_MR_result.tsv"), header = T, sep = "\t")
  print(filename)
  print(res)
  
  Visualizing_MR_forest(res,save_plot = TRUE,plot_pdf = "forest_plot")
}



