library(QTLMR)
library(readxl)
library(stringr)  

finn_dir = "/home/stardust/Documents/finngen_metabolism_mtag"
finn_files = list.files(path = finn_dir, pattern = ".txt", full.names = T)
metadata = read_excel("/home/stardust/Documents/finngene_summary_table.xlsx")
omopids = str_extract(basename(finn_files), "(?<=R12_)\\d+(?=_MTAG\\.txt)")
use_omopids = c("3009041", "3024561", "3026493", "40758310", "44786774")
selected_finn_files = finn_files[which(omopids %in% use_omopids)]
save_paths = paste0("/home/stardust/Documents/HDL_L_result/", use_omopids)
for (i in 2:length(selected_finn_files)) {
    if (!dir.exists(save_paths[i])) {dir.create(save_paths[i])}
    dat_rg <- HDL_rg_parallel(GWASfile = c("/home/stardust/Documents/msa_gwas_formatfile/msa_24_MTAG.txt",selected_finn_files[i]),
                            LD.path = "/home/stardust/Documents/HDL/UKB_imputed_SVD_eigen99_extraction",
                            Nref = 335265,
                            N0 = NULL,
                            numCores = 24,
                            eigen.cut = "automatic",
                            jackknife.df = FALSE,
                            intercept.output = FALSE,
                            fill.missing.N = NULL,
                            lim = exp(-18),
                            verbose = FALSE,
                            save_name = use_omopids[i],
                            save_path = save_paths[i])
}