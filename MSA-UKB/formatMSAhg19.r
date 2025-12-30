library(QTLMR)

# MSA GWAS PATHS:
msapath1 = "/home/stardust/Documents/MSA-UKB-GWAS/GCST90406924/GCST90406924.h_hg19.txt"
msapath2 = "/home/stardust/Documents/MSA-UKB-GWAS/GCST90406925/GCST90406925.h_hg19.txt"
# SNP	CHR	BP	effect_allele	other_allele	beta	standard_error	effect_allele_frequency	p_value	rsid	
# ci_upper	ci_lower	REF	n	hm_coordinate_conversion	odds_ratio	hm_code	variant_id

format_dat(
    msapath1,
    type = "outcome",
    snps=NULL,
    header = TRUE,
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "standard_error",
    eaf_col = "effect_allele_frequency",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "p_value",
    samplesize_col = "n",
    min_pval = 1e-200,
    chr_col = "CHR",
    pos_col = "BP",
    log_pval = FALSE,
    Twosample_dat = TRUE,
    SMR_dat = FALSE,
    MTAG_dat = FALSE,
    METAL_dat = FALSE,
    GWAS_name = "GCST90406924_hg19",
    save_path = "/home/stardust/Documents/MSA-UKB-GWAS/"
)

format_dat(
    msapath2,
    type = "outcome",
    snps=NULL,
    header = TRUE,
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "standard_error",
    eaf_col = "effect_allele_frequency",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "p_value",
    samplesize_col = "n",
    min_pval = 1e-200,
    chr_col = "CHR",
    pos_col = "BP",
    log_pval = FALSE,
    Twosample_dat = TRUE,
    SMR_dat = FALSE,
    MTAG_dat = FALSE,
    METAL_dat = FALSE,
    GWAS_name = "GCST90406925_hg19",
    save_path = "/home/stardust/Documents/MSA-UKB-GWAS/"
)