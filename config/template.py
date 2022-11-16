"""An area to assign project-specific variables. This is used alongside rhesus.yaml."""


NULL = ""

# Project directories to deposit output files
project = __file__.split("/")[-1].split(".")[0] + "/"  # Take project name from file name
location = "/master/username/variant-analysis/"  # Path to variant-analysis directory
resources = location + "resources/" + project
results = location + "results/" + project
log = location + "log/" + project

# Directory containing FASTQ files
reads = resources + "reads/"

# Only for GBS. A tab-delimited file with the first column as {organism_id} and the second as the barcode (e.g., AACGTT)
GBS_barcodes = resources + "barcodes/GBS_barcodes.tsv"

# Only for AMP. A tab-delimited file with the first column as {organism_id},
# the second as the i7 adapter, and the third as i5 adapter.
AMP_barcodes = resources + "barcodes/AMP_barcodes.tsv"

# One sample name per line
samples = resources + "samples/GBS_samples.list"

# FASTA reference genome
ref_fasta = resources + "ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz"  # Taken from ENSEMBL. Uses chromosomes as 1, 2 ... X, Y...

BQSR_known_variants = resources + "ref_vcf/macaca_mulatta.vcf.gz"


# These are for a single chromsome
SNP_VQSR_truth_vcf = resources + "truth_and_training/GbS_Common_snps79typed.recode.vcf.gz"
SNP_VQSR_truth_prior = 15.0
SNP_VQSR_training_vcf = resources + "truth_and_training/marmoset_gbs_pass_filter_vep_recal_beagle.vcf.gz"
SNP_VQSR_training_prior = 10.0
SNP_VQSR_known_vcf = NULL
SNP_VQSR_known_prior = 2.0
SNP_VQSR_truth_vcf = resources + "ref/mult_obs_Mmul_10.vcf.gz"
SNP_VQSR_truth_prior = 15.0
SNP_VQSR_training_vcf = resources + "ref/sing_obs_Mmul_10.vcf.gz"
SNP_VQSR_training_prior = 13.0
SNP_VQSR_known_vcf = NULL
SNP_VQSR_known_prior = 2.0

indel_VQSR_truth_vcf = NULL
indel_VQSR_truth_prior = 15.0
indel_VQSR_training_vcf = NULL
indel_VQSR_training_prior = 13.0
indel_VQSR_known_vcf = NULL
indel_VQSR_known_prior = 2.0
contig_remapping = NULL
indel_chromosome_remapping = NULL  # Inverse of contig_remapping

gene_counts = {
    "rna_transcripts": NULL,  # Just example. Doesn't exist.
    # A mapping that can be generated from a features_table file from NCBI
    "transcript_to_gene": NULL,  # Just example. Doesn't exist.
}

# Genotyped .vcf file
ref_vcf = NULL
# Tab-delimited file with 3 columns.
# First: Sample IDs used in .vcf.
# Second: Population IDs (for ancestry determination)
# Third: Individual IDs to be written into .geno file. May be same as first column.
pop_ids = NULL
# Number of prinicple components to use in PCA. I.e., dimension of reference space.
dim = 2

# For `relations.smk`
# Input .vcf file for lcMLkin
vcf = NULL
# List of founders for lcMLkin
founders = NULL
# Must have same name before .pop suffix as other PLINK input files to do supervised admixture
pop = NULL  # CURRENTLY DOESN'T USE THIS LINE, BUT MUST BE PLACED IN APPROPRIATE DIRECTORY. SEE CORRESPONDING RULE.
# Frequency below which to not consider as an AIM.
# Must be greater than 0.5, perferably as close to 0.9999 while still providing variants.
AIMs_max_diff = 0.85
# fst
gzvcf = NULL
# All populations for pairwise analysis
# Directory containing files with ids grouped by file.
pop_dir = NULL
# Lists of ids for each population. Includes those for Mmul_10.
pops = NULL

# PEDSYS-generated pedigree
pedigree = NULL
map = NULL
parents_table =None
sex_table = None
update_fids = None
effective_size: 8000

covergae = {
    # Minimum read depth to include for all samples.
    "cutoff": 0.9,
}

# Files to be generated (not currently using)
target_files = [
    "README.md",
    # User defined input goes below
    # These are two examples of how to structure desired files
    ## expand("{results}vcf/{sample}.g.vcf.gz", results=config["results"], sample=SAMPLE_NAMES)
    ## config["results"] + "db/cohort.sample_map"
    #expand(config["results"] + "genotypes/filtered/GBS_WES_WGS_labels.SNP.filtered.vcf.gz", chr=RHESUS_CHROMOSOMES),
    "alignments/raw/GBS16356.bam",
    ]