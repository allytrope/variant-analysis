# Create list of chromosomes. May require running Snakemake once to create these files and then again to use them in rules.
contigs_file = config["resources"] + "ref_fna/contigs.list"
if Path(contigs_file).is_file():
    with open(contigs_file) as f:
        CONTIGS = f.read().splitlines()
else:
    CONTIGS = []

chromosomes_file = config["resources"] + "ref_fna/chromosomes.list"
if Path(chromosomes_file).is_file():
    with open(chromosomes_file) as f:
        CHROMOSOMES = f.read().splitlines()
else:
    CHROMOSOMES = []

autosomes_file = config["resources"] + "ref_fna/autosomes.list"
if Path(autosomes_file).is_file():
    with open(autosomes_file) as f:
        AUTOSOMES = f.read().splitlines()
else:
    AUTOSOMES = []


## Output shortcuts
# A test
test = config["results"] + "test/test.txt"

# Create post-processed BAMs
bams = expand(config["results"] + "alignments/markdup/{collect}.bam",
    collect=collect_samples(fmt="{batch}/{seq}{indiv}_{library}_{flowcell_lane}"))

# Create gVCFs through GATK by library
gvcfs = expand(config["results"] + "gvcf/{collect}.chr{chr}.g.vcf.gz",
    collect=collect_samples(fmt="{batch}/{seq}{indiv}_{library}"),
    chr=AUTOSOMES)

# Create GenomicsDB datastore. Also for when updating the datastore with new samples
datastore = expand(config["results"] + "db/{dataset}/completed/{chr}.txt",
    dataset=config['dataset'],
    chr=AUTOSOMES)

# Create joint-called VCF through GATK
joint_called = expand(config["results"] + "joint_call/polyallelic/{dataset}.chr{chr}.vcf.gz",
    dataset=config['dataset'],
    chr=AUTOSOMES)

filtered_vcf = expand(config["results"] + "genotypes/filtered/{dataset}.SNP.chr{chr}.vcf.gz",
    dataset=config['dataset'],
    chr=AUTOSOMES)

pass_vcf = expand(config["results"] + "genotypes/pass/{dataset}.all.SNP.chr{chr}.bcf.csi",
    dataset=config['dataset'],
    chr=AUTOSOMES)

left_join_vcf_WES_WGS = expand(config["results"] + "genotypes/pass/{dataset1}+{dataset2}_left_join.all.SNP.chr{chr}.bcf",
    dataset1="WES4",
    dataset2="WGS3",
    chr=AUTOSOMES)

annotated = expand(config["results"] + "genotypes/annotated/bcsq/{dataset}.{subset}.SNP.chr{chr}.bcf",
    dataset=config['dataset'],
    subset=config['subset'],
    chr=AUTOSOMES)

Clinvar_tsv = expand(config["results"] + "liftover/rheMac10ToHg38/annotated/VEP/{dataset}.{subset}.SNP.chr{chr}.tsv",
    dataset=config['dataset'],
    subset=config['subset'],
    chr=AUTOSOMES)

# Format required by PMx
king_relatedness_PMx = expand(config["results"] + "kinship/KING/{dataset}.{subset}.SNP.autosomal.PMx.matrix",
    dataset=config['dataset'],
    subset=config['subset'])

# fROH - use in ROH/fROH_bcftools.ipynb
froh = expand(config["results"] + "roh/bcftools/{dataset}.{subset}.chr{chrom}.RG.roh",
    dataset=config['dataset'],
    subset=config['subset'],
    chrom=AUTOSOMES)

#GCTA
inbreeding = expand(config["results"] + "inbreeding/GCTA/pass/{dataset}.{subset}.ibc",
    dataset=config['dataset'],
    subset=config['subset'],
    )

# Heterozygosity
heterozygosity = config["results"] + 'heterozygosity/gvcf_counts.het',

# Unsupervised ADMIXTURE
unsupervised_admixture = expand(config["results"] + "admixture/ADMIXTURE/unsupervised/{dataset}.{subset}.SNP.autosomal.seed{seed}.{clusters}.Q",
    dataset=config['dataset'],
    subset=config['subset'],
    clusters=[2,3,4,5,6,7],
    seed=[899, 900, 901, 902, 903, 904, 905, 906, 907, 908]
    ),