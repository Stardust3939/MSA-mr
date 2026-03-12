import os
import math
import pandas as pd

filepath = "/home/stardust/Documents/LDSC_result"
files = os.listdir(filepath)
for file in files[0:5]:
    with open(os.path.join(filepath, file), 'r') as f:
        lines = f.readlines()
        data = []
        for line in lines[-5:-3]:
            data.append(line.strip().split())
        df = pd.DataFrame(data)
        if df[["p"]].iloc[0] <= 0.5:
            print(df)