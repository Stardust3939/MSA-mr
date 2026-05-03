library(QTLMR)
library(TwoSampleMR)

library(parallel)
library(doParallel)

data_source = "/home/stardust/Documents/eqtlMR/eqtldata"

files = list.files(path = data_source, pattern = ".gz", full.names = T)

cl <- makeCluster(8)
registerDoParallel(cl)

foreach (i = 1:length(files), .combine = 'c', .packages = c('QTLMR', "readxl","stringr","reticulate")) %dopar% {
    filename = basename(files[i])
    save_name = str_replace(filename, ".tsv.gz", "")
    format_dat(
    files[i],
    type = "exposure",
    snps=NULL,
    header = TRUE,
    snp_col = "MarkerID",
    beta_col = "BETA",
    se_col = "SE",
    phenotype_col = "feature",
    eaf_col = "AF_Allele2",
    effect_allele_col = "Allele2",
    other_allele_col = "Allele1",
    pval_col = "p.value",
    samplesize_col = "N",
    min_pval = 1e-200,
    chr_col = "CHR",
    pos_col = "POS",
    log_pval = FALSE,
    Twosample_dat = TRUE,
    SMR_dat = FALSE,
    MTAG_dat = FALSE,
    METAL_dat = FALSE,
    GWAS_name = save_name,
    save_path = "/home/stardust/Documents/eqtlMR/eqtlformat"
)
}