library(QTLMR)
library(foreach)
library(doParallel)
library(parallel)
library(readxl)  
library(stringr)

finn_dir = "/home/stardust/Documents/finngene_gz"

file_list = list.files(finn_dir, pattern = "*.gz$", full.names = T)

data_info_finngen = read_excel("/home/stardust/Documents/finngene_summary_table.xlsx")

out_dir = "/home/stardust/Documents/finngen_metabolism_rsid/"

dbsnp38dir = "/home/stardust/Documents/GCF_000001405.40.gz"
dbsnp19dir = "/home/stardust/Documents/GCF_000001405.25.gz"

filenames = basename(file_list)

omopids = str_extract(filenames, "(?<=R12_)\\d+(?=\\.gz)")

cl <- makeCluster(30)
registerDoParallel(cl)

foreach (i = 1:length(file_list), .combine = 'c', .packages = c('QTLMR', "readxl","stringr","reticulate")) %dopar% {
  filename = filenames[i]
  omopid = omopids[i]
  output_name = paste0("finngen_R12_", omopid, "_snp.txt")
  use_miniconda("r45")
  tran_SNP_python(GWAS_file = file_list[i],
                  dbsnp_file = dbsnp38dir,
                  build = 38,
                  chr_col = '"#chrom"',
                  pos_col = "pos",
                  effect_allele_col = "alt",
                  other_allele_col = "ref",
                  keep_snp = TRUE,
                  exact_map = FALSE,
                  save_name = output_name,
                  save_path = out_dir)
}
gc()
stopCluster(cl)