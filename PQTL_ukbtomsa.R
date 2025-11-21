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

pqtl_dat = data_cispQTL_ukb_ppp_exp(clump = 1)
pqtl_dat = data_cispQT_decode_exp(clump = 1)
pqtl_dat = data_cispQT_fenland_exp()
# count elements in exposure column:
exposures = table(pqtl_dat$exposure)
pqtl_dat <- remove_MHC_data(dat = pqtl_dat,
                            chr_col = "chr.exposure",
                            pos_col = "pos.exposure",
                            MHC_start = 28477797,
                            MHC_end = 33448354)


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

msaoutcome <- read_outcome_data(
  "D:\\data\\GCST90406924\\GCST90406924.tsv",
  snps = NULL,
  sep = "\t",
  #phenotype_col = "Phenotype",
  snp_col = "rs_id",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  samplesize_col = "n",
  #gene_col = "gene",
  #id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE,
  chr_col = "chromosome",
  pos_col = "base_pair_location"
)



write.table(msaoutcome, file = "D:\\data\\MSA_outcome_add.txt", sep = "\t", quote = FALSE, row.names = FALSE)

outcome_data <- modified_proxy_1000G(
  pqtl_dat$SNP,
  GWASfile = "D:\\data\\MSA_outcome_add.txt",
  proxies = T,
  rsq = 0.8,
  kb = 5000,
  nsnp = 5000,
  maf_threshold = 0.3,
  bfile_1000G = "D:\\data\\1kgv3\\EUR"
)

dat <- harmonise_data(exposure_dat = pqtl_dat, outcome_dat = outcome_data)

res <- xQTL_mr(dat, FDR_method = "fdr", PVE = TRUE)
# resuk is res with exposure = GUSB:
resuk <- res[res$exposure == "GUSB", ]
# modify resuk data, change all GUSB to GUSB_ukb_ppp:
resuk$exposure <- "GUSB_ukb_ppp"
resuk$outcome <- "MSA"
resde <- res[res$exposure == "GUSB", ]
resde$exposure <- "GUSB_decode"
resde$outcome <- "MSA"
resfin <- res[res$exposure == "GUSB", ]
resfin$exposure <- "GUSB_finngene"
resfin$outcome <- "MSA"
res <- rbind(resuk, resde, resfin)
xQTL_volcano_plot(
  res,
  breaks = c(0.1, 0.2, 0.3, 0.4, 0.5),
  scale_color = c("gray", "red", "green"),
  FDR_method = "fdr",
  save_plot = FALSE,
  pdf_name = "finn_volcano_plot",
  width = 6,
  height = 5,
  save_path = "D:\\"
)

res = read_excel("D:\\pic\\pic331.xlsx")

xQTL_forest(
  res,
  pvalsig = "fdr",
  ci_col = "#000000",
  ci_fill = "#000000",
  ci_lwd = 1,
  xlim = c(0, 8),
  ticks_at = c(0, 4, 8),
  arrange_OR = FALSE,
  arrange_expnm = FALSE,
  save_plot = TRUE,
  plot_pdf = "pic331",
  width = 8,
  height = 5,
  save_path = "D:\\pic"
)
dev.off()
