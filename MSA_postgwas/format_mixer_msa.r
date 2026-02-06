library(QTLMR)
library(foreach)
library(doParallel)
library(parallel)
library(readxl)  
library(stringr)

msa1 = "/home/stardust/Documents/MSA-UKB-GWAS/GCST90406924/GCST90406924.h.tsv"
msa2 = "/home/stardust/Documents/MSA-UKB-GWAS/GCST90406925/GCST90406925.h.tsv"

format_dat(
    msa1,
    type = "exposure",
    snps = NULL,
    header = TRUE,
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "standard_error",
    eaf_col = "effect_allele_frequency",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    chr_col = "chromosome",
    pos_col = "base_pair_location",
    samplesize_col = "n",
    min_pval = 1e-200,
    Twosample_dat = FALSE,
    SMR_dat = FALSE,
    MTAG_dat = TRUE,
    METAL_dat = FALSE,
    GWAS_name = "msa_24",
    save_path = "/home/stardust/Documents"
)

format_dat(
    msa2,
    type = "exposure",
    snps = NULL,
    header = TRUE,
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "standard_error",
    eaf_col = "effect_allele_frequency",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    chr_col = "chromosome",
    pos_col = "base_pair_location",
    samplesize_col = "n",
    min_pval = 1e-200,
    Twosample_dat = FALSE,
    SMR_dat = FALSE,
    MTAG_dat = TRUE,
    METAL_dat = FALSE,
    GWAS_name = "msa_25",
    save_path = "/home/stardust/Documents"
)

gc()

msa_mtag1 = "/home/stardust/Documents/msa_24_MTAG.txt"
msa_mtag2 = "/home/stardust/Documents/msa_25_MTAG.txt"

MiXeR_format_dat(GWASfile = msa_mtag1,
                 ncase = 8016,
                 ncontrol = 0,
                 opt_arguments_csv = NULL,
                 opt_arguments_qc = NULL,
                 exclude_ranges = "6:26000000-34000000",
                 save_name = "msa_mixer1",
                 save_path = "/home/stardust/Documents/")

MiXeR_format_dat(GWASfile = msa_mtag2,
                 ncase = 8016,
                 ncontrol = 0,
                 opt_arguments_csv = NULL,
                 opt_arguments_qc = NULL,
                 exclude_ranges = "6:26000000-34000000",
                 save_name = "msa_mixer2",
                 save_path = "/home/stardust/Documents/")