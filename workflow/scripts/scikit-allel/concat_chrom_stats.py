import pickle as pk

import pandas as pd

# Snakemake variables
concat_pickle = snakemake.output.pickle
pickles = snakemake.input.pickles

# Load data
dfs = []
for pickle in pickles:
    chrom = list(filter(lambda x: "chr" in x, pickle.split(".")))[0].strip("chr")
    with open(pickle, "rb") as f:
        df = pk.load(f)
        df["chrom"] = chrom  # Create new column detailing chromosome
        dfs.append(df)
df = pd.concat(dfs)

df.to_pickle(concat_pickle)