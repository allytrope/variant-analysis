"""Calculate runs of homozygosity for each chromosome of each individual."""

import allel
import numpy as np
import pandas as pd
import polars as pl

# # Variables taken from Snakemake rule
# vcf = snakemake.input.vcf
# contig_lengths = snakemake.input.contig_lengths
# roh_pickle = snakemake.output.roh_pickle
# froh_pickle = snakemake.output.froh_pickle
# chromosome = snakemake.wildcards.chr

# Variables taken from Snakemake rule
from sys import argv
vcf = argv[1]
contig_lengths = argv[2]
roh_pickle = argv[3]
froh_pickle = argv[4]
chromosome = argv[5]

# Read contig lengths from TSV into dictionary
lengths = dict()
with open(contig_lengths, "r") as f:
    print(contig_lengths)
    for row in f:
        print("After")
        print(row)
        row = row.strip().split()
        lengths[row[0]] = int(row[1])

# Create empty roh list for all samples
roh_list = []

# Create empty froh dictionary, which will contain proportion of genome in a ROH for each sample.
froh_list = []

# Pull genotype data from VCF
callset = allel.read_vcf(vcf, ['variants/CHROM', 'variants/POS', 'samples', 'calldata/GT'])

# Loop through each sample
for sample_idx, sample in enumerate(callset['samples']):
    for chrom in [chromosome]: # sorted(list(set(callset['variants/CHROM']))):
        chrom_sites = np.equal(callset['variants/CHROM'], [chrom])

        # Select genotypes by chromosome and sample
        genotypes = allel.GenotypeArray(callset['calldata/GT'][chrom_sites])
        genotypes = genotypes[:, sample_idx]

        positions = callset['variants/POS'][chrom_sites]


        # Create is_accessible list
        # Reading BED into a list the length of genome
        #file = "/master/abagwell/variant-analysis/results/rhesus/coverage/common_WES_0.5_loci.bed"
        file = "/master/abagwell/variant-analysis/resources/rhesus/annotations/exons.tsv"
        bed = pl.read_csv(file, has_header=False, separator="\t", new_columns=["chrom", "start", "end"], dtypes=[pl.String, pl.Int32, pl.Int32])
        
        # Exonic positions marked as accessible
        exons = bed.filter(
            pl.col("chrom") == chromosome
        ).with_columns(
            position = pl.int_ranges("start", "end")
        ).drop(
            "start", "end"
        ).explode("position").with_columns(
            is_accessible = pl.lit(True)
        )
        
        # Genomic positions marked as inaccessible
        genome = pl.arange(0, lengths[chromosome], eager=True).to_frame("position").with_columns(
            # Add chrom column
            chrom = pl.lit(chromosome)
        ).select(
            # Reorder columns
            "chrom", "position"
        ).with_columns(
            is_accessible = pl.lit(False)
        )
         # Replace matching genomic rows with exonic ones
        is_accessible = list(genome.join(exons, on="position", how="left").with_columns(
            is_accessible = pl.col("is_accessible_right").fill_null(False)
        ).drop("chrom_right", "is_accessible_right")["is_accessible"])


        roh_df, froh = allel.roh_poissonhmm(genotypes, positions, phet_roh=0.00000043, min_roh=1_000_000, contig_size=lengths[chrom])
        #roh_df, froh = allel.roh_poissonhmm(genotypes, positions, phet_roh=0.0001, min_roh=1_000_000, contig_size=lengths[chrom], is_accessible=is_accessible)
        #roh_df, froh = allel.roh_mhmm(genotypes, positions, phet_roh=0.0001, min_roh=1_000_000, contig_size=lengths[chrom])


        # Append to samples list
        roh_df.insert(0, 'sample', sample)
        roh_df.insert(1, 'chrom', chrom)
        roh_list.append(roh_df)

        # Add froh to dictionary
        froh_list.append((sample, chrom, froh))

# Store data as pickle files.
# These can then be opened into a pandas dataframe and used to make graphs with Seaborn
#pd.concat(roh_list, ignore_index=True).to_pickle(roh_pickle)
#pd.DataFrame(froh_list, columns=['sample', 'chrom', 'froh']).to_pickle(froh_pickle)

pd.concat(roh_list, ignore_index=True).to_csv(roh_pickle, sep="\t", index=False)
pd.DataFrame(froh_list, columns=['sample', 'chrom', 'froh']).to_csv(froh_pickle, sep="\t", index=False)
