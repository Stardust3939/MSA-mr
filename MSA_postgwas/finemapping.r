library(QTLMR)

res <- FM_summary(GWASfile = "xxx_MTAG.txt",
                  bfile = "G:/QTLMR_test/1kg.v3/EUR",
                  index_snps = c("rs11966200", "rs145365399", "rs9391847"),
                  index_snp_wind_kb = 500,
                  R2_threshold = 0.4,
                  snp_col = "SNP",
                  chr_col = "CHR",
                  pos_col = "BP",
                  pval_col = "pval",
                  save_name = "index",
                  save_path = "./FM")