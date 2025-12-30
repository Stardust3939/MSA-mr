library(QTLMR)

# load MSA GWAS tsv file:
msa_gwas <- read.table("/home/stardust/Documents/MSA-UKB-GWAS/GCST90406925/GCST90406925.h.tsv", header = TRUE, sep = "\t")

# rename columns to standard names:
colnames(msa_gwas)[which(colnames(msa_gwas) == "chromosome")] <- "CHR"
colnames(msa_gwas)[which(colnames(msa_gwas) == "base_pair_location")] <- "BP"
colnames(msa_gwas)[which(colnames(msa_gwas) == "rs_id")] <- "SNP"


# liftover from hg38 to hg19:
dat <- GWAS_liftover(msa_gwas,      #注意此处的GWAS不是路径，是已读入至R语言中的数据
                     build_from = "hg38",
                     build_to = "hg19",
                     chain_file = NULL,
                     chain_source = "ensembl",
                     imputation_ind = FALSE,
                     chrom_col = "CHR",     
                     start_col = "BP",
                     as_granges = FALSE,
                     style = "NCBI",
                     verbose = TRUE,
                     save_name = "GCST90406925.h_hg19",
                     save_path = "/home/stardust/Documents/MSA-UKB-GWAS/GCST90406925/")