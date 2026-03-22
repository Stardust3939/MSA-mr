import pandas
import os

source_dir = "/home/stardust/Documents/finn_v2"
dir_list = os.listdir(source_dir)
# sort dir_list in ascending order
dir_list.sort()

for dir_name in dir_list:
    dir_path = os.path.join(source_dir, dir_name)
    # check if dir_path is a directory:
    if not os.path.isdir(dir_path):
        continue
    # list all files in dir_path:
    file_list = os.listdir(dir_path)
    # find file name in file_list ending with "MR_result.tsv":
    file_name = ""
    for f in file_list:
        if f.endswith("MR_result.tsv"):
            file_name = os.path.join(dir_path, f)
            break
    if not file_name:
        continue
    df = pandas.read_csv(file_name, sep="\t", low_memory=False)
    if 'pval' not in df.columns:
        continue
    # count how many rows in column 'pval' less than 0.05:
    count = len(df[df['pval'] < 0.05])
    if count >= 2:
        # find file name end with "heterogeneity.tsv" and "pleiotropy.tsv"
        heterogeneity_file = ""
        pleiotropy_file = ""
        for f in file_list:
            if f.endswith("heterogeneity.tsv"):
                heterogeneity_file = os.path.join(dir_path, f)
            if f.endswith("pleiotropy.tsv"):
                pleiotropy_file = os.path.join(dir_path, f)
        if heterogeneity_file and pleiotropy_file:
            heterogeneity_df = pandas.read_csv(heterogeneity_file, sep="\t", low_memory=False)
            pleiotropy_df = pandas.read_csv(pleiotropy_file, sep="\t", low_memory=False)
            # check if there is any row in heterogeneity_df with pval less than 0.05 and any row in pleiotropy_df with pval less than 0.05
            if (heterogeneity_df['Q_pval'] < 0.05).any() and (pleiotropy_df['pval'] < 0.05).any():
                print(dir_name, count, "heterogeneity pval < 0.05 and pleiotropy pval < 0.05")
            elif (heterogeneity_df['Q_pval'] < 0.05).any():
                print(dir_name, count, "heterogeneity pval < 0.05")
            elif (pleiotropy_df['pval'] < 0.05).any():
                print(dir_name, count, "pleiotropy pval < 0.05")
            else:
                print(dir_name, count)
