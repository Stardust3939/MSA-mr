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

for (i in 1:length(selected_finn_files)) {
    MiXeR_run(GWASfile = c(msa_files[1], selected_finn_files[i]),
          GWASname = c("MSA", use_omopids[i]),
          bfile = "/home/stardust/Documents/mixer1000g",
          mixer_py = "/home/stardust/Documents/gsa-mixer/precimed/mixer.py",
          mixer_figures_py = "/home/stardust/Documents/gsa-mixer/precimed/mixer_figures.py",
          lib = "/home/stardust/Documents/gsa-mixer/src/build/lib/libbgmg.so",
          statistic = "mean std",
          fit1_diffevo_fast_repeats = 20,
          fit2_diffevo_fast_repeats = 20,
          save_path = save_paths[i],
          cores = 30)
}