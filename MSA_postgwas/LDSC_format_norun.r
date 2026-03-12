#############################################
#############################################
## This file is only a record of the code  ##
## It is not intended to be run as a whole ##
#############################################
#############################################

library(QTLMR)
library(foreach)
library(doParallel)
library(parallel)
library(readxl)  
library(stringr)

# first format finngen gz data to TwoSampleMR format:
# read metadata
finn_dir = "/home/stardust/Documents/finngen_metabolism_rsid_n"
finn_files = list.files(path = finn_dir, pattern = ".txt", full.names = T)
finn_name = basename(finn_files)

# format finngen data to TwoSampleMR format
cl <- makeCluster(8)
registerDoParallel(cl)

foreach (i = 1:length(finn_files), .combine = 'c', .packages = c('QTLMR', "readxl","stringr","reticulate")) %dopar% {
    # if destination file already exists, skip
    if (file.exists(paste0("/home/stardust/Documents/finn_tmr/",finn_name[i], "_TwosampleMR.txt"))) {
        print(paste0("File ", finn_name[i], " already exists, skipping..."))
        next
    }
    format_dat(
    finn_files[i],
    type = "exposure",
    snps=NULL,
    header = TRUE,
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "sebeta",
    eaf_col = "af_alt",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    pval_col = "pval",
    samplesize_col = "N",
    min_pval = 1e-200,
    chr_col = "#chrom",
    pos_col = "pos",
    log_pval = FALSE,
    Twosample_dat = TRUE,
    SMR_dat = FALSE,
    MTAG_dat = FALSE,
    METAL_dat = FALSE,
    GWAS_name = finn_name[i],
    save_path = "/home/stardust/Documents/finn_tmr"
)
    gc()
}
gc()
stopCluster(cl)

format_dat(
    "/home/stardust/Documents/MSA-UKB-GWAS/GCST90406925/GCST90406925.h.tsv",
    type = "exposure",
    snps=NULL,
    header = TRUE,
    snp_col = "rs_id",
    beta_col = "beta",
    se_col = "standard_error",
    eaf_col = "effect_allele_frequency",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "p_value",
    samplesize_col = "n",
    min_pval = 1e-200,
    chr_col = "chromosome",
    pos_col = "base_pair_location",
    log_pval = FALSE,
    Twosample_dat = TRUE,
    SMR_dat = FALSE,
    MTAG_dat = FALSE,
    METAL_dat = FALSE,
    GWAS_name = "msa25",
    save_path = "/home/stardust/Documents"
)

LDSC_py_sumstats(GWASfile = "/home/stardust/Documents/msa_gwas_formatfile/msa25_TwosampleMR.txt",
                GWAS_name = "msa25",
                N = 8016,
                MAF_min = 0.01,
                opt_arguments = NULL,
                merge_alleles = "/home/stardust/Documents/LDSC/w_hm3.noMHC.snplist",
                save_path = "/home/stardust/Documents")

finn_dir = "/home/stardust/Documents/finn_tmr"
finn_files = list.files(path = finn_dir, pattern = ".txt", full.names = T)
metadata = read_excel("/home/stardust/Documents/finngene_summary_table.xlsx")
# extract omopids from file names: *R12_<OMOPID>_snp*
omopids = str_extract(basename(finn_files), "(?<=R12_)\\d+(?=_snp)")
use_omopids = c("3009041", "3024561", "3026493", "40758310", "44786774")
selected_finn_files = finn_files[which(omopids %in% use_omopids)]
save_paths = paste0("/home/stardust/Documents/LDSC_result/", use_omopids)



for (i in 1:length(selected_finn_files)) {
    LDSC_py_sumstats(GWASfile = selected_finn_files[i],
                GWAS_name = use_omopids[i],
                N = metadata$num_cases[which(metadata$OMOPID == use_omopids[i])],
                MAF_min = 0.01,
                opt_arguments = NULL,
                merge_alleles = "/home/stardust/Documents/LDSC/w_hm3.noMHC.snplist",
                save_path = "/home/stardust/Documents/LDSC_result/")
}