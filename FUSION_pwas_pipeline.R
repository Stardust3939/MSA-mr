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
#library(cause)
library(vroom)
library(MungeSumstats)
#library(GenomicFiles)
#library(meta)
#library(readr)
#library(readxl)
#library(forestploter)
library(ldscr)
setwd("D:/data/FUSION_PWAS")
TWAS_fusion_format_data(
  GWASfile = "D:\\data\\MSA_outcome_add.txt",
  N = NULL,
  save_name = "MSA_add",
  save_path = "./"
)

TWAS_fusion_assoc_test(
  test_help = FALSE,
  Sumstatsfile = "./MSA_add.sumstats",
  weights = "./PWAS_EA/Plasma_Protein_EA_hg38.pos",
  weights_dir = "./PWAS_EA/Plasma_Protein_weights_EA",
  resource = "other",
  ref_ld_chr = "./LDREF/1000G.EUR.",
  start_chr = 1,
  end_chr = 22,
  coloc_pval = NULL,
  GWASN = NULL,
  perm = 10000,
  PANELN = NA,
  opt_arguments = NULL,
  remove_MHC = TRUE,
  FDR_method = "bonferroni",
  save_name = "MSAadd",
  save_path = "./MSAadd",
  cores = 12
)

TWAS_fusion_assoc_test(
  test_help = FALSE,
  Sumstatsfile = "./MSA.sumstats",
  weights = "./ROSMAP.n376.fusion.WEIGHTS/train_weights.pos.addN",
  weights_dir = "./ROSMAP.n376.fusion.WEIGHTS",
  resource = "other",
  ref_ld_chr = "./LDREF/1000G.EUR.",
  start_chr = 1,
  end_chr = 22,
  coloc_pval = NULL,
  GWASN = NULL,
  perm = 10000,
  PANELN = NA,
  opt_arguments = NULL,
  remove_MHC = TRUE,
  FDR_method = "bonferroni",
  save_name = "rosmap",
  save_path = "./rosmap",
  cores = 12
)
TWAS_fusion_assoc_test(
  test_help = FALSE,
  Sumstatsfile = "./MSA.sumstats",
  weights = "./Banner.n152.fusion.WEIGHTS/train_weights.pos",
  weights_dir = "./Banner.n152.fusion.WEIGHTS",
  resource = "other",
  ref_ld_chr = "./LDREF/1000G.EUR.",
  start_chr = 1,
  end_chr = 22,
  coloc_pval = NULL,
  GWASN = NULL,
  perm = 10000,
  PANELN = NA,
  opt_arguments = NULL,
  remove_MHC = TRUE,
  FDR_method = "bonferroni",
  save_name = "banner",
  save_path = "./banner",
  cores = 12
)
TWAS_fusion_assoc_test(
  test_help = FALSE,
  Sumstatsfile = "./MSA.sumstats",
  weights = "./dir02_PWASweights/TWAS_Weights.pos",
  weights_dir = "./dir02_PWASweights",
  resource = "other",
  ref_ld_chr = "./LDREF/1000G.EUR.",
  start_chr = 1,
  end_chr = 22,
  coloc_pval = NULL,
  GWASN = NULL,
  perm = 10000,
  PANELN = NA,
  opt_arguments = NULL,
  remove_MHC = TRUE,
  FDR_method = "bonferroni",
  save_name = "banner",
  save_path = "./banner",
  cores = 12
)

mdata = read.table("D:\\data\\FUSION_PWAS\\SCZ\\SCZ.csv", header = TRUE, sep = ",")

# keep ID,CHR column in mdata:
mdata1 = mdata[,c("ID","CHR","BEST.GWAS.ID","TWAS.P")]

write.table(mdata1, "D:\\data\\FUSION_PWAS\\SCZ\\SCZ2.csv", sep = ",", row.names = FALSE, col.names = TRUE)

tran_chrpos_from_SNP(
  GWASfile = "D:\\data\\FUSION_PWAS\\SCZ\\SCZ2.csv",
  snp_col = "BEST.GWAS.ID",
  chr_col = "chr.outcome",
  pos_col = "pos.outcome",
  build = 38,
  save_path = "D:\\data\\FUSION_PWAS\\SCZ",
  save_name = "mdata2"
)
mdata2 = read.table("D:\\data\\FUSION_PWAS\\SCZ\\mdata2.txt", header = TRUE, sep = "\t")
# change mdata2 column sequence:
mdata2 = mdata2[,c("SNP","CHR.x","pos.outcome","TWAS.P")]
Visualizing_Manhattan(mdata2,
                      plot.type="m",
                      LOG10=TRUE,
                      ylim=NULL,
                      col = "grey",
                      threshold=c(1e-8,1e-4),
                      threshold.lty=c(1,2),
                      threshold.lwd=c(1,1),
                      threshold.col=c("black","grey"),
                      amplify=FALSE,
                      bin.size=1e6,
                      chr.den.col=c("darkgreen", "yellow", "red"),
                      signal.col=c("red","green"),
                      signal.cex=c(1.5,1.5),
                      signal.pch=c(19,19),
                      highlight="rs9530",
                      highlight.text="GUSB",
                      highlight.text.font = 1,
                      file="jpg",
                      file.name="test",
                      dpi=300,
                      file.output=TRUE,
                      verbose=TRUE,
                      width=14,
                      height=6)



