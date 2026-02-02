library(QTLMR)
library(foreach)
library(doParallel)
library(parallel)
library(readxl)  
library(stringr)

gwas_dir = "/home/stardust/Documents/finngen_metabolism_rsid"
file_list = list.files(gwas_dir, pattern = "*.txt$", full.names = T)
file_names = basename(file_list)
omopids = str_extract(file_names, "(?<=R12_)\\d+(?=_snp\\.txt)")

cl <- makeCluster(30)
registerDoParallel(cl)

foreach (i = 1:length(file_list), .combine = 'c', .packages = c('QTLMR', "readxl","stringr")) %dopar% {
  filename = file_names[i]
  omopid = omopids[i]
  output_name = paste0("finngen_R12_", omopid)
  # if output file already exists, skip
  output_file = paste0("/home/stardust/Documents/finngen_metabolism_mtag/", output_name, "_MTAG.txt")
  if (file.exists(output_file)) {
    cat(paste0("File ", output_file, " already exists. Skipping...\n"))
    next
  }
  format_dat(
    file_list[i],
    type = "exposure",
    snps = NULL,
    header = TRUE,
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "sebeta",
    eaf_col = "af_alt",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    pval_col = "pval",
    chr_col = "#chrom",
    pos_col = "pos",
    min_pval = 1e-200,
    Twosample_dat = FALSE,
    SMR_dat = FALSE,
    MTAG_dat = TRUE,
    METAL_dat = FALSE,
    GWAS_name = output_name,
    save_path = "/home/stardust/Documents/finngen_metabolism_mtag/"
)
gc()
}
stopCluster(cl)

