library(QTLMR)
library(TwoSampleMR)

library(parallel)
library(doParallel)

datasource = "/home/stardust/Documents/eqtlMR/eqtlformat"

files = list.files(path = datasource, pattern = ".txt", full.names = T)


cl <- makeCluster(16)
registerDoParallel(cl)

foreach (i = 1:length(files), .combine = 'c', .packages = c('QTLMR', "readxl","stringr","reticulate","TwoSampleMR")) %dopar% {
    filename = basename(files[i])
    save_name = str_replace(filename, ".txt", "")
    dat <- read.table(files[i], header = TRUE, sep = "\t")
    dat <- remove_MHC_data(dat,
                chr_col = "chr.exposure",
                pos_col = "pos.exposure",
                MHC_start = 28000000,
                MHC_end = 34000000)
    
    clump_dat <- clump_data_local_Online(dat,
                               snp_col = "SNP",
                               pval_col = "pval.exposure",
                               clump_method = "PVAL",
                               beta_col = "beta.exposure",
                               se_col = "se.exposure",
                               clump_pval = 5e-05,
                               clump_kb = 10000,
                               clump_r2 = 0.001,
                               pop = "EUR",
                               bfile_1000G = "/home/stardust/Documents/1kgv3/EUR")
    
    rm(dat)
    gc()
    outcome_dat <- modified_proxy_1000G(snps = clump_dat$SNP,
                                    GWASfile = "/home/stardust/Documents/msa_gwas_formatfile/msa24_outcome_TwosampleMR.txt",
                                    proxies = T,
                                    rsq = 0.8,
                                    kb = 5000,
                                    nsnp = 5000,
                                    maf_threshold = 0.01,
                                    bfile_1000G = "/home/stardust/Documents/1kgv3/EUR")

    h_dat <- harmonise_data(exposure_dat = clump_dat, outcome_dat = outcome_dat, action = 2)
    rm(clump_dat, outcome_dat)
    gc()
    f_dat <- calculation_Fvalue_R2(exp_GWAS = h_dat, samplesize = NULL)
    rm(h_dat)
    gc()
    # f_dat contains Fvalue_Scheme1 to Fvalue_Scheme4, get max Fvalue among the four schemes for each SNP
    f_dat$Fvalue <- apply(f_dat[, c("Fvalue_Scheme1", "Fvalue_Scheme2", "Fvalue_Scheme3", "Fvalue_Scheme4")], 1, max, na.rm = TRUE)
    # filter out SNPs with Fvalue < 10
    f_dat <- f_dat[f_dat$Fvalue >= 10, ]
    res<- xQTL_mr(f_dat, FDR_method = "fdr", PVE = TRUE)
    rm(f_dat)
    gc()
    res <- xQTL_volcano_plot(res,
                         breaks = c(0.1, 0.2, 0.3, 0.4, 0.5),
                         scale_color = c("green", "grey", "red"),
                         FDR_method = "fdr",
                         save_plot = TRUE,
                         pdf_name = save_name,
                         width = 6,
                         height = 5,
                         save_path = "/home/stardust/Documents/eqtlMR/results")
    forestname = paste0(save_name, "_forest")
    xQTL_forest(res,
            pvalsig = "fdr",
            ci_col = "#000000",
            ci_fill = "#000000",
            ci_lwd = 1,
            xlim = c(0, 6),
            ticks_at = c(0, 0.5, 1, 2, 3, 4, 5, 6),
            save_plot = TRUE,
            plot_pdf = forestname,
            width = 8,
            height = 5,
            save_path = "/home/stardust/Documents/eqtlMR/results")
}