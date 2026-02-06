library(QTLMR)
library(readxl)
library(stringr)  
Sys.setenv(R_MAX_MEM_SIZE = "90000M")
finn_dir = "/home/stardust/Documents/finngen_metabolism_mixer"
finn_files = list.files(path = finn_dir, pattern = ".gz", full.names = T)
msa_dir = "/home/stardust/Documents/msa_gwas_formatfile"
msa_files = list.files(path = msa_dir, pattern = ".gz", full.names = T)

metadata = read_excel("/home/stardust/Documents/finngene_summary_table.xlsx")
omopids = str_extract(basename(finn_files), "(?<=R12_)\\d+(?=_qc_noMHC\\.csv\\.gz)")

# use omopids:3009041 3024561 3026493 40758310 44786774
use_omopids = c("3009041", "3024561", "3026493", "40758310", "44786774")
selected_finn_files = finn_files[which(omopids %in% use_omopids)]

save_paths = paste0("/home/stardust/Documents/mixer_results_selected/", use_omopids)

for (i in 2:length(selected_finn_files)) {
    MiXeR_Univariate (GWASfile = selected_finn_files[i],
                  bfile = "/home/stardust/Documents/mixer1000g",
                  mixer_py = "/home/stardust/Documents/gsa-mixer/precimed/mixer.py",
                  lib = "/home/stardust/Documents/gsa-mixer/src/build/lib/libbgmg.so",
                  stars_rep = 1,
                  end_rep = 20,
                  fit1 = TRUE,
                  opt_arguments_fit = NULL,
                  opt_arguments_test = NULL,
                  save_path = save_paths[i],
                  cores = 30)
}
