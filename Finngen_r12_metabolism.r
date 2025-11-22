library(QTLMR)
library(TwoSampleMR)
df = load("Z:\\finngen_R12_1e_4.Rdata")

exposure_list = unique(data$exposure)

msaoutcome = read_exposure_data(
  "Z:\\MSA_TwosampleMR.txt",
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
list = read_excel("Z:\\finngene_summary_table.xlsx")

setwd("Z:\\Finngen_r12_metabolism_MR_results")


for (finn_exposure in exposure_list) {
  exposure_dat = data[data$exposure == finn_exposure, ]
  colnames(exposure_dat)[colnames(exposure_dat) == "SNP"] <- "rsid"
  colnames(exposure_dat)[colnames(exposure_dat) == "pval.exposure"] <- "pval"
  colnames(exposure_dat)[colnames(exposure_dat) == "id.exposure"] <- "trait_id"

  exposure_clump <- ieugwasr::ld_clump(dat = exposure_dat, clump_r2 = 0.01, clump_kb = 10000, clump_p = 1, plink_bin = "C:/Users/Jerry/.conda/envs/r45/Lib/R/library/plinkbinr/bin/plink_Windows.exe", bfile = "Z:/1kgv3/EUR")

  colnames(exposure_clump)[colnames(exposure_clump) == "rsid"] <- "SNP"
  colnames(exposure_clump)[colnames(exposure_clump) == "pval"] <- "pval.exposure"
  colnames(exposure_clump)[colnames(exposure_clump) == "trait_id"] <- "id.exposure"

  # get exposure number from list, finn_exposure is LONGNAME, get OMOPID:
  exposure_index = which(list$LONGNAME == finn_exposure)
  filename = paste0(list$OMOPID[exposure_index])
  # create output folder for each exposure:
  if (!dir.exists(filename)) {
    dir.create(filename)
  }
  setwd(paste0("Z:\\Finngen_r12_metabolism_MR_results\\", filename))
  

  dat <- harmonise_data(exposure_dat = exposure_clump, outcome_dat = msaoutcome)
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