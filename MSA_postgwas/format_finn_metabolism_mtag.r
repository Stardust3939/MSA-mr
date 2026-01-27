library(QTLMR)
library(foreach)
library(doParallel)
library(parallel)
library(readxl)  
library(stringr)

finn_dir = "/mnt/nas_ssd/finngene_gz"

file_list = list.files(finn_dir, pattern = "*.gz$", full.names = T)

data_info_finngen = read_excel("/mnt/nas_ssd/finngene_summary_table.xlsx")

out_dir = "/mnt/nas_ssd/finngen_metabolism_mtag/"

filenames = basename(file_list)

omopids = str_extract(filenames, "(?<=R12_)\\d+(?=\\.gz)")
file_path = file_list[1]
mopid = omopids[1]
out_path = paste0(out_dir, omopid)
dat_finngen <- format_data_FinnGen(GWASfile= file_path,
                                   GWAS_name = omopid,
                                   save_path = out_path,
                                   type = "outcome",
                                   min_pval = 1e-200,
                                   build_to_hg19 = TRUE,
                                   Twosample_dat = TRUE,
                                   SMR_dat = TRUE,
                                   MTAG_dat = TRUE,
                                   METAL_dat = TRUE)

