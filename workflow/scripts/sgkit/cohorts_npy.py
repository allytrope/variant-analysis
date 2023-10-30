import numpy as np
import zarr

# Load populations into a list of sets
pops = []
for file in snakemake.input.pops:
    with open(file) as f:
        pops.append(set(f.read().strip().split("\n")))

# Load zarr
ds = zarr.load(snakemake.input.zarr)

# Create list of integers to tell which population each sample belongs to
# Same sample order as in VCF
cohorts = []
for sample in ds["sample_id"]:
    assigned_flag = False
    for index, pop in enumerate(pops):
        if sample in pop:
            if assigned_flag == True:
                raise Exception(f"Sample '{sample}' assigned to multiple populations.")
            cohorts.append(index)
            assigned_flag = True
    if assigned_flag == False:
        cohorts.append(-1)
np.save(snakemake.output.npy, np.array(cohorts, dtype=np.int16))