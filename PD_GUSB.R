

setwd("D:/data/PD_GUSB")


format_data_IEU_VCF(
  IEU_VCF = "./ieu-b-7.vcf.gz",
  type = "outcome",
  ncase = 33674,
  ncontrol = 449056,
  samplesize = 482730,
  min_pval = 1e-200,
  Twosample_dat = TRUE,
  SMR_dat = FALSE,
  MTAG_dat = FALSE,
  METAL_dat = FALSE,
  save_name = "PD",
  save_path = "./"
)

outcome_data = read.table("./PD_TwosampleMR.txt", header = T, sep = "\t")

pqtl_dat1 = data_cispQTL_ukb_ppp_exp(clump = 1)
pqtl_dat2 = data_cispQT_decode_exp(clump = 1)
pqtl_dat3 = data_cispQT_fenland_exp()

pqtl_dat1 <- remove_MHC_data(dat = pqtl_dat1,
                             chr_col = "chr.exposure",
                             pos_col = "pos.exposure",
                             MHC_start = 28477797,
                             MHC_end = 33448354)
pqtl_dat2 <- remove_MHC_data(dat = pqtl_dat2,
                             chr_col = "chr.exposure",
                             pos_col = "pos.exposure",
                             MHC_start = 28477797,
                             MHC_end = 33448354)

pqtl_dat3 <- remove_MHC_data(dat = pqtl_dat3,
                             chr_col = "chr.exposure",
                             pos_col = "pos.exposure",
                             MHC_start = 28477797,
                             MHC_end = 33448354)

pqtl_dat1 = pqtl_dat1[pqtl_dat1$exposure == "GUSB",]
pqtl_dat2 = pqtl_dat2[pqtl_dat2$exposure == "GUSB",]
pqtl_dat3 = pqtl_dat3[pqtl_dat3$exposure == "GUSB",]


file = file3
filename = basename(file)
pqtl_dat = pqtl_dat2

outcome_data <- modified_proxy_1000G(
  pqtl_dat$SNP,
  GWASfile = "./PD_TwosampleMR.txt",
  proxies = T,
  rsq = 0.8,
  kb = 10000,
  nsnp = 5000,
  maf_threshold = 0.1,
  bfile_1000G = "D:\\data\\1kgv3\\EUR"
)

dat <- harmonise_data(exposure_dat = pqtl_dat, outcome_dat = outcome_data)

res2 <- xQTL_mr(dat, FDR_method = "fdr", PVE = TRUE)

res1$exposure <- "GUSB_UKB"
res2$exposure <- "GUSB_decode"
res <- rbind(res1, res2)

xQTL_volcano_plot(
  res1,
  breaks = c(0.1, 0.2, 0.3, 0.4, 0.5),
  scale_color = c("green", "grey", "red"),
  FDR_method = "fdr",
  save_plot = FALSE,
  pdf_name = "volcano_plot",
  width = 6,
  height = 5,
  save_path = "D:\\"
)

xQTL_forest(
  res,
  pvalsig = "fdr",
  ci_col = "#000000",
  ci_fill = "#000000",
  ci_lwd = 1,
  xlim = c(0.8, 1.2),
  ticks_at = c(0.8, 0.9, 1, 1.1, 1.2),
  arrange_OR = FALSE,
  arrange_expnm = FALSE,
  save_plot = FALSE,
  plot_pdf = "forest_plot",
  width = 12,
  height = 5,
  save_path = "D:\\"
)
# delete nsnp column in res:
res$nsnp <- NULL
