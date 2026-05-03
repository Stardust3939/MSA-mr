library(QTLMR)

res_Multi <- TWAS_fusion_Multi_test(Sumstatsfile = "/home/stardust/Documents/msa_gwas_formatfile/msa24.sumstats.gz",  
                                    weights_type = "nofilter.pos",  
                                    weights_dir = "/home/stardust/Documents/Fusion_TWAS/WEIGHTS",  
                                    resource = "GTEx_v8", 
                                    ref_ld_chr = "/home/stardust/Documents/Fusion_TWAS/LDREF/1000G.EUR.", 
                                    start_chr = 1, 
                                    end_chr = 22, 
                                    coloc_pval = 0.05,  
                                    GWASN = 8016, 
                                    perm = 10000,
                                    PANELN = NA, 
                                    opt_arguments = NULL, 
                                    remove_MHC = TRUE, 
                                    FDR_method = "fdr",  
                                    save_name = "MSA24", 
                                    save_path = "/home/stardust/Documents/Fusion_TWAS/results", 
                                    cores = 24)