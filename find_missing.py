import pandas
import os

dir = "Z:\\Finngen_r12_metabolism_MR_results"

list = os.listdir(dir)
filelist = pandas.read_excel("Z:\\finngene_summary_table.xlsx")

filelist = filelist['OMOPID']

# find missing id in list
for name in list:
    if name not in filelist:
        print(name)

list = pandas.DataFrame(list, columns=["OMOPID"])

list.to_csv("Z:\\Finngenlist.csv", index=False)