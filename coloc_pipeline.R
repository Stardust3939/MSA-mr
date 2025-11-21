setwd("D:\\")
iinfo = data_info_ukb_ppp_pQTL()
# find GUSB in info:
gusb_info = iinfo[iinfo$HGNC.symbol == "GUSB", ]

gusb_data <- get_cispQTL_decode_ukb_Online(
  gene = "GUSB",
  gene_ID = "OID20296",
  resource = "ukb_ppp",
  cis_wind_kb = 1000,
  build = 38,
  save_name = "GUSB",
  save_path = "./",
  SMR_data = TRUE
)

decode_info = data_info_decode_pQTL()
head(decode_info)
gusb_info = decode_info[decode_info$gene_name_Ensembl == "GUSB", ]
gusb_info
gusb_data <- get_cispQTL_decode_ukb_Online(
  gene = "GUSB",
  gene_ID = "15562_24",
  resource = "decode",
  cis_wind_kb = 1000,
  build = 38,
  save_name = "GUSB_decode",
  save_path = "./",
  SMR_data = TRUE
)

fenland_info = data_info_fenland_pQTL()
head(fenland_info)
gusb_info = fenland_info[fenland_info$gene_name == "GUSB", ]
gusb_info

gusb_data <- get_cispQTL_fenland_Online(
  gene = "GUSB",
  SomaScan_ids = "15562_24",
  cis_wind_kb = 1000,
  build = 38,
  save_name = "GUSB_fenland",
  save_path = "./",
  SMR_data = TRUE
)

# add MAF column to pqtl_dat, if eaf.exposure<= 0.5, MAF = eaf.exposure, else MAF = 1 - eaf.exposure:
pqtl_dat$MAF <- ifelse(pqtl_dat$eaf.exposure <= 0.5, pqtl_dat$eaf.exposure, 1 - pqtl_dat$eaf.exposure)

# rename pqtl_dat columns, origin name: chr.exposure; new name: CHR:
colnames(pqtl_dat)[colnames(pqtl_dat) == "chr.exposure"] <- "CHR"
# pval.exposure to pval:
colnames(pqtl_dat)[colnames(pqtl_dat) == "pval.exposure"] <- "pval"
# samplesize.exposure to N:
colnames(pqtl_dat)[colnames(pqtl_dat) == "samplesize.exposure"] <- "N"
# beta.exposure to beta:
colnames(pqtl_dat)[colnames(pqtl_dat) == "beta.exposure"] <- "beta"
# se.exposure to se:
colnames(pqtl_dat)[colnames(pqtl_dat) == "se.exposure"] <- "se"
# pos.exposure to BP:
colnames(pqtl_dat)[colnames(pqtl_dat) == "pos.exposure"] <- "BP"

# check if pqtl_dat contains duplicates in SNP, if contains, remove duplicates:
if (any(duplicated(pqtl_dat$SNP))) {
  pqtl_dat <- pqtl_dat[!duplicated(pqtl_dat$SNP), ]
}

#keep only chr = 7 in pqtl_dat:
pqtl <- pqtl_dat[pqtl_dat$CHR == 7, ]
# keep only BP between 64000000 and 67000000 in pqtl:
pqtl <- pqtl[pqtl$BP >= 60000000 & pqtl$BP <= 70000000, ]

pqtl = read.table("D:\\cis_GUSB_15562_24.txt", header = TRUE, sep = "\t")
#pqtl = pqtl_dat[2708:2711,]
coloc_GWAS(
  pqtl,
  msaoutcome,
  type_exposure = "quant",
  col_pvalues_exposure = "pval",
  col_N_exposure = "N",
  col_MAF_exposure = "MAF",
  col_beta_exposure = "beta",
  col_se_exposure = "se",
  col_snp_exposure = "SNP",
  col_chr_exposure = "CHR",
  col_pos_exposure = "BP",
  sd_exposure = NA,
  type_outcome = "cc",
  col_pvalues_outcome = "pval.outcome",
  col_N_outcome = "samplesize.outcome",
  col_MAF_outcome = NA,
  col_beta_outcome = "beta.outcome",
  col_se_outcome = "se.outcome",
  col_snp_outcome = "SNP",
  prevalence_outcome = NA,
  save_stacked_dat = FALSE,
  build = 38,
  save_locus = TRUE,
  title1 = "xQTL",
  title2 = "GWAS",
  width = 8,
  height = 5,
  plot_pdf = "locuscompare_fenland",
  save_path = "./"
)

#check if rs34356500 in msaoutcome$SNP, output the entire row if rs34356500 in msaoutcome$SNP:
if ("rs34356500" %in% msaoutcome$SNP) {
  print("rs34356500 in msaoutcome")
  print(msaoutcome[msaoutcome$SNP == "rs34356500", ])
} else {
  print("rs34356500 not in msaoutcome")
}
if ("rs34356500" %in% outcome_data$SNP) {
  print("rs34356500 in outcome_data")
} else {
  print("rs34356500 not in outcome_data")
}


