"""An area to assign project-specific variables."""

from pathlib import Path

from snakemake.io import expand


## These may be left unchanged.
NULL = ""
project = __file__.split("/")[-1].split(".")[0] + "/"  # Take project name from file name

## Users should add files and paths here.
# Project directories to deposit output files
path = "/master/username/variant-analysis/"  # Path to variant-analysis directory
resources = path + "resources/" + project
results = path + "results/" + project
log = path + "log/" + project

# Directory containing FASTQ files
reads = resources + "reads/"

barcodes = {
    # Only for GBS. A tab-delimited file with the first column as {organism_id} and the second as the barcode (e.g., AACGTT)
    "GBS_barcodes": resources + "barcodes/GBS_barcodes.tsv",
    
    # Only for AMP. A tab-delimited file with the first column as {organism_id},
    # the second as the i7 adapter, and the third as i5 adapter.
    "AMP_barcodes": resources + "barcodes/AMP_barcodes.tsv",
}

# One sample name per line
samples = resources + "samples/GBS_samples.list"

# FASTA reference genome
ref_fasta = resources + "ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz"  # Taken from ENSEMBL. Uses chromosomes as 1, 2 ... X, Y...

BQSR_known_variants = resources + "ref_vcf/macaca_mulatta.vcf.gz"

# Tab-delimited. Column1 is child id, column2 is sire id, column3 is dam id. Also, remove any names with underscores.
pedigree = resources + "pedigree/pedigree.tsv",

filtering = {
    "min_AF": 0.05,
    "max_AF": 0.95,
    "min_AC": 10,
}

## Variables to use in target files (and used in .smk files)
# Create list of samples
with open(samples) as f:
    SAMPLE_NAMES = f.read().splitlines()
print(SAMPLE_NAMES)

# Create list of chromosomes
if Path(results + "db/chromosomes.list").is_file():
    with open(results + "db/chromosomes.list") as f:
        CHROMOSOMES = f.read().splitlines()
    print(CHROMOSOMES)
else:
    print("No chromosome file yet.")

# Files to create. The paths can be found as under "output" in the rules of .smk files.
# Every output should start with "results" to direct to the correct path.
# A simple target file:
# results + "haplotypes/SHAPEIT4/all_animals.SNP.chr17.vcf.gz",
# Example of a more complex target file:
# expand(results + "vcf/{sample}.g.vcf.gz", sample=SAMPLE_NAMES)
target_files = [
]
