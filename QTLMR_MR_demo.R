library(QTLMR)
library(TwoSampleMR)

format_dat(
  msadata,
  type = "outcome",
  snps = NULL,
  header = TRUE,
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect.allele.exposure",
  other_allele_col = "other.allele.exposure",
  pval_col = "pval_origin.exposure",
  chr_col = "chr.exposure",
  Twosample_dat = TRUE,
  GWAS_name = "MSA",
  save_path = "D:\\"
)

pqtl_dat = data_cispQTL_ukb_ppp_exp(clump = 2)

pqtl_dat <- remove_MHC_data(dat = pqtl_dat,
                           chr_col = "chr.exposure",
                           pos_col = "pos.exposure",
                           MHC_start = 28477797,
                           MHC_end = 33448354)

# 整理为结局
dat <- format_dat(
  ukb,
  type = "outcome",
  snps = NULL,
  header = TRUE,
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  id_col = "id.exposure",
  gene_col = "gene.exposure",
  chr_col = "chr.exposure",
  info_col = "info.exposure",
  min_pval = 1e-200,
  log_pval = FALSE,
  Twosample_dat = TRUE,
  SMR_dat = FALSE,
  MTAG_dat = FALSE,
  METAL_dat = FALSE,
  GWAS_name = "UKB",
  save_path = "D:\\"
)

# 读取UKB-ppp结局数据
ukb_outcome = read.table("D:\\UKB_TwosampleMR.txt", header = TRUE, sep = "\t")

dat <- harmonise_data(exposure_dat = msadata, outcome_dat = ukb_outcome)


dat <- harmonise_data(exp_dat,outcome_dat)


#计算F值与R2

dat <- calculation_Fvalue_R2(dat)

#进行MR分析

res <- xQTL_mr(dat, FDR_method = "fdr", PVE = TRUE)

res <- xQTL_volcano_plot(res,
                         breaks = c(0.1, 0.2, 0.3, 0.4, 0.5),
                         scale_color = c("green", "grey", "red"),
                         FDR_method = "fdr",
                         save_plot = TRUE,
                         pdf_name = "volcano_plot",
                         width = 6,
                         height = 5,
                         save_path = "./outdata/")
#绘制森林图
xQTL_forest(res,
            pvalsig = "fdr",
            ci_col = "#000000",
            ci_fill = "#000000",
            ci_lwd = 1,
            xlim = c(0, 3.5),
            ticks_at = c(0, 0.5, 1, 2, 3, 3.5),
            save_plot = TRUE,
            plot_pdf = "forest_plot",
            width = 8,
            height = 5,
            save_path = "./outdata/")

# 获取ukb_ppp数据
