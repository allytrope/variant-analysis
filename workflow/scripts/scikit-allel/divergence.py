import allel
import numpy as np
import pandas as pd

pops = snakemake.input.pops
pickle = snakemake.output.pickle
vcf = snakemake.input.vcf  #"/master/abagwell/variant-analysis/results/rhesus/haplotypes/SHAPEIT5_WGS/SNPRC_WGS_WES.SNP.chr18.vcf.gz"
window_size = snakemake.wildcards.window_size  #50_000

callset = allel.read_vcf(vcf, ['variants/POS', 'samples', 'calldata/GT'])

def indices_of_pop_members(population, callset):
    """Map the population members to indices of sample ids in the VCF."""
    with open(population, "r") as f:
        pop = f.read().strip().split("\n")

    sample_indices = []
    samples = list(callset["samples"])
    for sample in pop:
        try:
            idx = samples.index(sample)
            sample_indices.append(idx)
        except:
            print(f"{sample} not in VCF")
    return sorted(sample_indices)

pop1 = indices_of_pop_members(pops[0], callset)
pop2 = indices_of_pop_members(pops[1], callset)

print(len(pop1))
print(len(pop2))

# Check for exclusivity of populations
intersection = set(pop1).intersection(set(pop2))
if intersection:
    raise Exception(f"The following samples are found in both populations: {intersection}")

pops = {"pop1": pop1, "pop2": pop2}

# Take genotypes from current chromosome
g = allel.GenotypeArray(callset['calldata/GT'])

# Allele counts
ac = g.count_alleles_subpops(pops)

# Calculate divergence
dxy, windows, n_bases, counts = allel.windowed_divergence(callset['variants/POS'], ac["pop1"], ac["pop2"], size=int(window_size)*1000, start=1, step=int(window_size)*1000//2)

# Create DataFrame
df = pd.DataFrame({
    'dxy': dxy,
    'start': np.hsplit(windows, 2)[0].ravel(),
    'stop': np.hsplit(windows, 2)[1].ravel(),
    'n_bases': n_bases,
    'counts': counts
})

# Pickle
df.to_pickle(pickle)