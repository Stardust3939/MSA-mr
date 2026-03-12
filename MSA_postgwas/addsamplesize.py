import pandas
import re
import os
from multiprocessing import Process

source = "/home/stardust/Documents/finngen_metabolism_rsid"
dest = "/home/stardust/Documents/finngen_metabolism_rsid_n"

filenames = os.listdir(source)
metadata ="/home/stardust/Documents/finngene_summary_table.xlsx"


def worker(metadatadest,sourcefold,destfold,filenames,startnum,tonum):
    metadata = pandas.read_excel(metadatadest)
    for i in range(startnum, tonum):
        filename = filenames[i]
        omopid = re.search(r"(?<=R12_)\d+(?=_snp\.txt\.txt)", filename).group(0)
        sample_size_row = metadata[metadata['OMOPID'] == int(omopid)]
        if not sample_size_row.empty:
            sample_size = sample_size_row['num_cases'].values[0]
            df = pandas.read_csv(os.path.join(sourcefold, filename), sep="\t")
            df['N'] = sample_size
            outname = "withsamplesize_" + filename
            df.to_csv(os.path.join(destfold, outname), sep="\t")

if __name__ == "__main__":
    num_processes = 30
    total_files = len(filenames)
    files_per_process = total_files // num_processes + 1
    processes = []

    for i in range(num_processes):
        startnum = i * files_per_process
        tonum = (i + 1) * files_per_process if i != num_processes - 1 else total_files
        p = Process(target=worker, args=(metadata, source, dest, filenames, startnum, tonum))
        processes.append(p)
        p.start()

    for p in processes:
        p.join()