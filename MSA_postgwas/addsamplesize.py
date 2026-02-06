import pandas
import re
import os

source = "/home/stardust/Documents/finngen_metabolism_mtag"
dest = "/home/stardust/Documents/finngen_metabolism_mtag_withsamplesize"

filenames = os.listdir(source)
metadata = pandas.read_excel("/home/stardust/Documents/finngene_summary_table.xlsx")

for filename in filenames:
    #filename like:finngen_R12_1175426_MTAG.txt
    omopid = re.search(r"(?<=R12_)\d+(?=_MTAG\.txt)", filename).group(0)
    sample_size_row = metadata[metadata['OMOPID'] == int(omopid)]
    if not sample_size_row.empty:
        sample_size = sample_size_row['num_cases'].values[0]
        '''
        with open(os.path.join(source, filename), 'r') as infile, open(os.path.join(dest, filename), 'w') as outfile:
            header = infile.readline().strip()
            outfile.write(header + '\n')
            for line in infile:
                parts = line.strip().split(',')
                parts.append(str(sample_size))
                outfile.write(','.join(parts) + '\n')
                '''
        print(f"Processed {filename} with sample size {sample_size}")