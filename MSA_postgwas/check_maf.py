import pandas
import re

source = "/home/stardust/Documents/uaieu.ma"
# read as tsv files:
df = pandas.read_csv(source, sep="\t")
# check if column "MAF" contains values not between 0 and 1
maf_values = df["MAF"]
invalid_maf = maf_values[(maf_values < 0) | (maf_values > 1)]
if not invalid_maf.empty:
    print("Invalid MAF values found:")
    print(invalid_maf)

# check if p values are between 0 and 1
p_values = df["P"]
invalid_p = p_values[(p_values < 0) | (p_values > 1)]
if not invalid_p.empty:
    print("Invalid P values found:")
    print(invalid_p)

# check if beta contains 0
beta_values = df["BETA"]
zero_beta = beta_values[beta_values == 0]
if not zero_beta.empty:
    print("Beta values equal to 0 found:")
    print(zero_beta)