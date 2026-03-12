library(QTLMR)
library(foreach)
library(doParallel)
library(parallel)
library(readxl)  
library(stringr)

finn_dir = "/home/stardust/Documents/finn_sumstat"
finn_files = list.files(path = finn_dir, pattern = ".gz", full.names = T)

omopids = str_extract(basename(finn_files), "\\d+(?=\\.sumstats\\.gz)")




for (i in 1:length(finn_files)) {
    LDSC_py_rg(test_help = FALSE,
           Sumstatsfile = c("/home/stardust/Documents/msa_gwas_formatfile/msa24.sumstats.gz",
                            finn_files[i]),
           ref_ld_chr = "/home/stardust/Documents/LDSC/eur_w_ld_chr",
           w_ld_chr = "/home/stardust/Documents/LDSC/eur_w_ld_chr",
           opt_arguments = NULL,
           save_name = omopids[i],
           save_path = "/home/stardust/Documents/LDSC_result/")
}




