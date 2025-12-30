library(QTLMR)
library(TwoSampleMR)
library(ieugwasr)
library(foreach)
library(doParallel)
library(parallel)
library(readxl)  

msapath1 = "/home/stardust/Documents/MSA-UKB-GWAS/GCST90406924_hg19_TwosampleMR.txt"
msapath2 = "/home/stardust/Documents/MSA-UKB-GWAS/GCST90406925_hg19_TwosampleMR.txt"
result_path = "/home/stardust/Documents/MSA-UKB-GWAS/result1/"
exposure_p_val_threshold = 5e-8
clump_r2_threshold = 0.001
clump_kb_threshold = 10000

msaoutcome <- read_outcome_data(
  msapath1,
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

files <- list.files("/home/stardust/Documents/MSA-UKB-GWAS/MR-GWAS-format", pattern = ".txt", full.names = T)

cl <- makeCluster(16)
registerDoParallel(cl)

foreach (i = 1:28, .combine = 'c', .packages = c('QTLMR', "TwoSampleMR","ieugwasr","readxl")) %dopar% {
  filename = basename(files[i])
  dir.create(paste0(result_path, filename, "_MR"))
  setwd(paste0(result_path, filename, "_MR"))

  # chr.exposure	pos.exposure	SNP	other_allele.exposure	effect_allele.exposure	beta.exposure	se.exposure	pval.exposure	
  # eaf.exposure	samplesize.exposure	id.exposure	exposure	mr_keep.exposure	pval_origin.exposure
  
  exposure_dat <- read_exposure_data(
    filename = files[i],
    clump = FALSE,
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
  
  exposure_dat_pval <- exposure_dat[exposure_dat$pval.exposure <= exposure_p_val_threshold,]
  
  #rename columns SNP to rsid in dataframe msa:
  colnames(exposure_dat_pval)[colnames(exposure_dat_pval) == "SNP"] <- "rsid"
  colnames(exposure_dat_pval)[colnames(exposure_dat_pval) == "pval.exposure"] <- "pval"
  colnames(exposure_dat_pval)[colnames(exposure_dat_pval) == "id.exposure"] <- "trait_id"
  
  
 # msa_clump <- clump_data(msa, clump_r2 = 0.01, clump_kb = 5000, pop = "EUR")
  
  
  
  exposure_clump <- ieugwasr::ld_clump(dat = exposure_dat_pval, clump_r2 = clump_r2_threshold, clump_kb = clump_kb_threshold, clump_p = 1, plink_bin = "/home/stardust/miniconda3/envs/r45/bin/plink", bfile = "/home/stardust/Documents/MSA-UKB-GWAS/1kgv3/EUR")
  
  colnames(exposure_clump)[colnames(exposure_clump) == "rsid"] <- "SNP"
  colnames(exposure_clump)[colnames(exposure_clump) == "pval"] <- "pval.exposure"
  colnames(exposure_clump)[colnames(exposure_clump) == "trait_id"] <- "id.exposure"
  
  
  dat <- harmonise_data(exposure_dat = exposure_clump, outcome_dat = msaoutcome)
  
  res <- mr(dat)
  #export res to tsv in working directory:
  wd = getwd()
  write.table(res, file = paste0(wd, "/", filename, "_MR_result.tsv"), sep = "\t", row.names = F)
  
  het <- mr_heterogeneity(dat)
  #export het to tsv in working directory:
  write.table(het, file = paste0(wd, "/", filename, "_MR_heterogeneity.tsv"), sep = "\t", row.names = F)
  
  ple <- mr_pleiotropy_test(dat)
  #export ple to tsv in working directory:
  write.table(ple, file = paste0(wd, "/", filename, "_MR_pleiotropy.tsv"), sep = "\t", row.names = F)
  
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

stopCluster(cl)
