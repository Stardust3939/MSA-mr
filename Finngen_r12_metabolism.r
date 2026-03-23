library(QTLMR)
library(TwoSampleMR)
library(readxl)
library(foreach)
library(doParallel)
library(parallel)
library(stringr)  



msaoutcome = read_outcome_data(
  "/home/stardust/Documents/msa_gwas_formatfile/msa25_TwosampleMR.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  min_pval = 1e-200,
  chr_col = "chr.exposure",
  pos_col = "pos.exposure",
  log_pval = FALSE
)

# open detailed exposure list to get exposure index:
# list = read_excel("Z:\\finngene_summary_table.xlsx")

# setwd("Z:\\Finngen_r12_metabolism_MR_results")
exposure_dir = "/mnt/nas_ssd/finn_tmr"

exposure_files = list.files(exposure_dir, pattern = "*.txt", full.names = TRUE)

cl <- makeCluster(8)
registerDoParallel(cl)

foreach (i = 1:383, .combine = 'c', .packages = c('QTLMR', "TwoSampleMR", "readxl", "stringr")) %dopar% {
  finn_exposure = exposure_files[i]
  filename = basename(finn_exposure)
  omopid = str_extract(filename, "(?<=R12_)\\d+(?=_snp)")

  outpath = paste0("/home/stardust/Documents/finn_25v2/", omopid)
  # if outpath does not exist, create it,else skip to next iteration:
  if (!dir.exists(outpath)) {
    dir.create(outpath)
  } else {
    next
  }
  setwd(outpath)

  exposure_dat = read_exposure_data(
    filename = finn_exposure,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure",
    pval_col = "pval.exposure",
    samplesize_col = "samplesize.exposure",
    min_pval = 1e-200,
    chr_col = "chr.exposure",
    pos_col = "pos.exposure",
    log_pval = FALSE
  )
  # filter exposure_dat by pval.exposure < 5e-6:
  exposure_dat_filter = subset(exposure_dat, pval.exposure < 5e-5)
  N = exposure_dat_filter[1,"samplesize.exposure"]
  exposure_dat_filter=transform(exposure_dat_filter,R2=2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))
  exposure_dat_filter=transform(exposure_dat_filter,F=(N-2)*R2/(1-R2))
  #filter exposure_dat by F > 10:
  exposure_dat_filter = subset(exposure_dat_filter, F > 10)

  ndata = nrow(exposure_dat_filter)
  if (ndata < 2) {
    next
  }
  # write an empty file to the outpath to indicate that this omopid has been processed:
  file.create(paste0(outpath, "/F_processed.txt"))
  rm(exposure_dat)
  gc()

  colnames(exposure_dat_filter)[colnames(exposure_dat_filter) == "SNP"] <- "rsid"
  colnames(exposure_dat_filter)[colnames(exposure_dat_filter) == "pval.exposure"] <- "pval"
  colnames(exposure_dat_filter)[colnames(exposure_dat_filter) == "id.exposure"] <- "trait_id"

  exposure_clump <- ieugwasr::ld_clump_local(dat = exposure_dat_filter, clump_r2 = 0.001, clump_kb = 10000, clump_p = 1, plink_bin = "/home/stardust/miniconda3/envs/r45/bin/plink", bfile = "/mnt/nas_ssd/1kgv3/EUR")
  gc()
  colnames(exposure_clump)[colnames(exposure_clump) == "rsid"] <- "SNP"
  colnames(exposure_clump)[colnames(exposure_clump) == "pval"] <- "pval.exposure"
  colnames(exposure_clump)[colnames(exposure_clump) == "trait_id"] <- "id.exposure"

  sizeofexposure = nrow(exposure_clump)
  if (sizeofexposure < 2) {
    next
  }
  # write an empty file to the outpath to indicate that this omopid has been clumped:
  file.create(paste0(outpath, "/clumped.txt"))
  rm(exposure_dat_filter)
  gc()

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
  
  rm(dat,exposure_dat,exposure_dat_filter,exposure_clump,res,het,ple,result_single,result_loo,p1,p2,p3,p4,)
  gc()
}

stopCluster(cl)
rm(list = ls())
gc()