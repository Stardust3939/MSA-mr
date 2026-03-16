import os
import math
import pandas as pd
from io import StringIO

filepath = "/home/stardust/Documents/LDSC_result"
files = os.listdir(filepath)
for file in files:
    with open(os.path.join(filepath, file), 'r') as f:
        lines = f.readlines()
        data = []
        for line in lines[-5:-3]:
            data.append(line.strip().split())
        df = pd.read_csv(StringIO('\n'.join(['\t'.join(row) for row in data])), sep='\t')
        if df['p'][0] < 0.05:
            # select column rg, se, z and p
            print(os.path.basename(file))
            print(df[['rg', 'se', 'z', 'p']])