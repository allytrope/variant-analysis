import allel
import pandas as pd

# Snakemake variables
pickle = snakemake.output.pickle
subpop = snakemake.input.subpop
vcf = snakemake.input.vcf

callset = allel.read_vcf(vcf, ['variants/POS', 'samples', 'calldata/GT'])

# Read samples file
with open(subpop, "r") as f:
    subpop_samples = f.read().strip().split("\n")

samples = list(callset["samples"])

# Find indices of samples in VCF
sample_indices = []
for sample in subpop_samples:
    try:
        sample_indices.append(samples.index(sample))
    except:
        pass
sample_indices = sorted(sample_indices)

# Subset data
subpop_data = callset["calldata/GT"].take(sample_indices, axis=1)

g = allel.GenotypeArray(subpop_data)
af = g.count_alleles().to_frequencies()

expected = allel.heterozygosity_expected(af, ploidy=2)

#np.save(output.npy, expected, allow_pickle=False, fix_imports=False)

# Create DataFrame
df = pd.DataFrame({
    'position': callset['variants/POS'],
    'exp_heterozygosity': expected,
})

# Pickle
df.to_pickle(pickle)