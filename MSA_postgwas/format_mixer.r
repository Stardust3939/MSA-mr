library(QTLMR)
library(foreach)
library(doParallel)
library(parallel)
library(readxl)  
library(stringr)

mtag_dir = "/home/stardust/Documents/finngen_metabolism_mtag_withsamplesize"
file_list = list.files(mtag_dir, pattern = "*.txt$", full.names = T)
file_names = basename(file_list)
omopids = str_extract(file_names, "(?<=R12_)\\d+(?=_MTAG\\.txt)")
cl <- makeCluster(20)
registerDoParallel(cl)

foreach (i = 3:length(file_list), .combine = 'c', .packages = c('QTLMR', "readxl","stringr")) %dopar% {
  filename = file_names[i]
  omopid = omopids[i]
  output_name = paste0("finngen_R12_", omopid)

  # open metadata excel file and find ncases for this omopid:
  metadata = read_excel("/home/stardust/Documents/finngene_summary_table.xlsx")
  ncases = metadata$`num_cases`[metadata$OMOPID == as.numeric(omopid)]
  
  # note that MiXeR_format_dat() WILL NOT fill in ncase or ncontrol if GWASfile["N"] is EMPTY!!!
  # FILL IN ncase AND ncontrol BEFORE running format!!! 
  MiXeR_format_dat(GWASfile = file_list[i],
                 ncase = ncases,
                 ncontrol = 0,
                 opt_arguments_csv = NULL,
                 opt_arguments_qc = NULL,
                 exclude_ranges = "6:26000000-34000000",
                 save_name = output_name,
                 save_path = "/home/stardust/Documents/finngen_metabolism_mixer/")
gc()
}
stopCluster(cl)