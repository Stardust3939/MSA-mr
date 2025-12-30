library(QTLMR)
library(ieugwasr)
library(readr)

folderpath = "/home/stardust/Documents/MSA-UKB-GWAS/MR-GWAS/"
outpath = "/home/stardust/Documents/MSA-UKB-GWAS/MR-GWAS-format/"
files = list.files(folderpath, pattern = "*.vcf.gz", full.names = TRUE)

# get gwas id, gwas id is the file name without suffix
gwas_ids = sub(".vcf.gz$", "", basename(files))
# print(gwas_ids)

# metadata_list = gwasinfo(id = gwas_ids, opengwas_jwt = get_opengwas_jwt())
# write_tsv(metadata_list, file = "/home/stardust/Documents/MSA-UKB-GWAS/GWAS_info.tsv")

# read meta.tsv in folderpath for gwas metadata
metadata_list = read_tsv("/home/stardust/Documents/MSA-UKB-GWAS/GWAS_info.tsv")
# go through each gwas, read sample_size from metadata_list:

for (i in 1:length(gwas_ids)) {
  gwas_id = gwas_ids[i]
  sample_size = metadata_list$sample_size[metadata_list$id == gwas_id]
  # full file name:
  file = files[i]
  dat_VCF <- format_data_IEU_VCF(IEU_VCF = file,
                               type = "exposure",
                               ncase = NA,
                               ncontrol = NA,
                               samplesize = sample_size,
                               min_pval = 1e-200,
                               Twosample_dat = TRUE,
                               SMR_dat = FALSE,
                               MTAG_dat = FALSE,
                               METAL_dat = FALSE,
                               save_name = gwas_id,
                               save_path = outpath)
}