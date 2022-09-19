import pandas as pd
from matminer.datasets import load_dataset


df = load_dataset("expt_gap_kingsbury")

print(df)

df2 = df.dropna()
print(df2)


df2.to_csv("matminer_mp.csv",sep=",")
