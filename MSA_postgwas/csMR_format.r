library(QTLMR)
library(csMR)

# record of command lines
# not intended for running, but for record keeping

csMR_format_dat(gwas_file = "/home/stardust/Documents/msa_gwas_formatfile/msa24_TwosampleMR.txt", 
                save_name = "msa24", 
                save_path = "/home/stardust/Documents/msa_gwas_formatfile/")

# serum UA
csMR_format_dat(gwas_file = "/mnt/data/finngen_formats/finn_tmr/withsamplesize_finngen_R12_3026493_snp.txt.txt_TwosampleMR.txt", 
                save_name = "uricacid", 
                save_path = "/home/stardust/Documents/")

csMR_format_dat(gwas_file = "/home/stardust/Documents/Uric_acid_19_TwosampleMR.txt", 
                save_name = "uaieu", 
                save_path = "/home/stardust/Documents/")

csMR_run(exposure_gwas = list(list(id = "UA",
                                   path = "/home/stardust/Documents/uaieu.ma",
                                   type = "quant")),
         sceqtl_input = list(list(id = "OG",
                                  path = "/home/stardust/Documents/csMR/sceQTL/Oligodendrocyte.ma")),
         outcome_dir = "/home/stardust/Documents/msa_gwas_formatfile/outcomes",
         reference_path = "/home/stardust/Documents/csMR/reference_genome_1000G_EUR",
         gwas_duplicated_snp_path = "None",
         coloc_window_size_bp = 100000,
         coloc_coverages = 0.9,
         coloc_threads = 30,
         coloc_cutoff = 0.8,
         cores = 30,
         save_path = "/home/stardust/Documents/csMR/results")