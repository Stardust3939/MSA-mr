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

for (finn_exposure in exposure_list) {
  exposure_dat = data[data$exposure == finn_exposure, ]
  colnames(exposure_dat)[colnames(exposure_dat) == "SNP"] <- "rsid"
  colnames(exposure_dat)[colnames(exposure_dat) == "pval.exposure"] <- "pval"
  colnames(exposure_dat)[colnames(exposure_dat) == "id.exposure"] <- "trait_id"

  exposure_clump <- ieugwasr::ld_clump(dat = exposure_dat, clump_r2 = 0.01, clump_kb = 10000, clump_p = 1, plink_bin = "C:/Users/Jerry/.conda/envs/r45/Lib/R/library/plinkbinr/bin/plink_Windows.exe", bfile = "Z:/1kgv3/EUR")

  colnames(exposure_clump)[colnames(exposure_clump) == "rsid"] <- "SNP"
  colnames(exposure_clump)[colnames(exposure_clump) == "pval"] <- "pval.exposure"
  colnames(exposure_clump)[colnames(exposure_clump) == "trait_id"] <- "id.exposure"
  
  harmonised_data <- harmonise_data(
    exposure_dat,
    msaoutcome,
    action = 2
  )
  
  mr_results <- mr(harmonised_data)
  
  output_file <- paste0("Z:\\Finngen_r12_metabolism_results\\", finn_exposure, "_MR_results.txt")
  write.table(mr_results, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
}