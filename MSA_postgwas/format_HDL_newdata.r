library(QTLMR)
library(foreach)
library(doParallel)
library(parallel)
library(readxl)  
library(stringr)

dat_VCF <- format_data_IEU_VCF(IEU_VCF = "/home/stardust/Documents/ebi-a-GCST90018977.vcf.gz",
                               type = "exposure",
                               ncase = NA,
                               ncontrol = NA,
                               samplesize = 343836,
                               min_pval = 1e-200,
                               Twosample_dat = FALSE,
                               SMR_dat = FALSE,
                               MTAG_dat = TRUE,
                               METAL_dat = FALSE,
                               save_name = "Uric_acid_19",
                               save_path = "/home/stardust/Documents")
# rename columns chr and pos
colnames(dat_VCF)[which(colnames(dat_VCF) == "chr.exposure")] <- "CHR"
colnames(dat_VCF)[which(colnames(dat_VCF) == "pos.exposure")] <- "BP"

dat <- GWAS_liftover(dat_VCF,
                     build_from = "hg19",
                     build_to = "hg38",
                     chain_file = "/home/stardust/Documents/GRCh37_to_GRCh38.chain",
                     chain_source = "ensembl",
                     imputation_ind = FALSE,
                     chrom_col = "CHR",     
                     start_col = "BP",
                     as_granges = FALSE,
                     style = "NCBI",
                     verbose = TRUE,
                     save_name = "Uric_acid_hg38",
                     save_path = "/home/stardust/Documents")

format_dat(
    dat,
    type = "exposure",
    snps=NULL,
    header = TRUE,
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    eaf_col = "eaf.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    pval_col = "pval.exposure",
    samplesize_col = "samplesize.exposure",
    min_pval = 1e-200,
    chr_col = "CHR",
    pos_col = "BP",
    log_pval = FALSE,
    Twosample_dat = FALSE,
    SMR_dat = FALSE,
    MTAG_dat = TRUE,
    METAL_dat = FALSE,
    GWAS_name = "uric_acid_hg38",
    save_path = "/home/stardust/Documents/"
)