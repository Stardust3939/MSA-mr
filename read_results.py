import pandas
import os

source_dir = "Z:\\Finngen_r12_metabolism_MR_results"
dir_list = os.listdir(source_dir)
# sort dir_list in ascending order
dir_list.sort()

for dir_name in dir_list:
    dir_path = os.path.join(source_dir, dir_name)
    file_name = os.path.join(dir_path, dir_name + "_MR_result.tsv")
    df = pandas.read_csv(file_name, sep="\t", low_memory=False)
    # count how many rows in column 'pval' less than 0.05:
    count = len(df[df['pval'] < 0.05])
    if count >= 3:
        print(dir_name, count)

