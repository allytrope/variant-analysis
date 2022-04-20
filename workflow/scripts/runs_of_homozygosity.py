"""Calculate runs of homozygosity for each chromosome of each individual."""

import allel
import numpy as np
import pandas as pd

# Variables taken from Snakemake rule
vcf = snakemake.input.vcf
roh_pickle = snakemake.output.roh_pickle
froh_pickle = snakemake.output.froh_pickle

# Reference genome Mmul_8.0.1 chromosomes
chromosome_lengths = {"1": 225_584_828,
               "2": 204_787_373,
               "3": 185_818_997,
               "4": 172_585_720,
               "5": 190_429_646,
               "6": 180_051_392,
               "7": 169_600_520,
               "8": 144_306_982,
               "9": 129_882_849,
               "10": 92_844_088,
               "11": 133_663_169,
               "12": 125_506_784,
               "13": 108_979_918,
               "14": 127_894_412,
               "15": 111_343_173,
               "16": 77_216_781,
               "17": 95_684_472,
               "18": 70_235_451,
               "19": 53_671_032,
               "20": 74_971_481,
               "X": 149_150_640}

# Create empty roh list for all samples
roh_list = []

# Create empty froh dictionary, which will contain proportion of genome in a ROH for each sample.
froh_list = []

# Pull genotype data from VCF
callset = allel.read_vcf(vcf, ['variants/CHROM', 'variants/POS', 'samples', 'calldata/GT'])

# Loop through each sample
for sample_idx, sample in enumerate(callset['samples']):
    for chrom in sorted(list(set(callset['variants/CHROM']))):
        chrom_sites = np.equal(callset['variants/CHROM'], chrom)

        # Select genotypes by chromosome and sample
        genotypes = allel.GenotypeArray(callset['calldata/GT'][chrom_sites])
        genotypes = genotypes[:, sample_idx]

        positions = callset['variants/POS'][chrom_sites]

        roh_df, froh = allel.roh_poissonhmm(genotypes, positions, min_roh=1_000_000, contig_size=chromosome_lengths[chrom])

        # Append to samples list
        roh_list.append(roh_df)

    # Add froh to dictionary
    froh_list.append((sample, chrom, froh))


# Store data as pickle files.
# These can then be opened into a pandas dataframe and used to make graphs with Seaborn
pd.concat(roh_list, ignore_index=True).to_pickle(roh_pickle)
pd.DataFrame(froh_list, columns=['samples', 'chrom', 'froh']).to_pickle(froh_pickle)

