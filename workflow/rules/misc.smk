rule create_readme:
    """Create README.md file in appropriate directory for easy reference about generated directories."""
    output: config["results"] + "README.md"
    shell: """echo '''
    Information about where these files fall in the workflow can be seen under the `.smk` files under `workflow/rules`.

    The directories for variant calling from `variant_calling.smk` and `phasing.smk` are created in the following order:

    ├─ alignments  # Duplicate sequences are marked and sorted
    |   |
    |   ├─ raw  # .bam files after .fastq files aligned to reference
    |   |
    |   ├─ markdup  # .bam files after running fixing mates pairs and marking duplicate reads
    |   |
    |   └─ alignments_recalibrated  # BQSR tables applied to sequences
    |       |
    |       └─ recal_tables  # BQSR recalibration tables that are applied to sequences
    |
    ├─ gvcf  # Called variants as .g.vcf files
    |
    ├─ db  # GenomicsDB datastore consolidating .g.vcf files. Not human-readable
    |
    ├─ joint_call  # db data converted to a joint .vcf file with VQSR applied
    |
    ├─ vcf  # Multisample VCF
    |   |
    |   ├─ hard_filtered # (Option 1) Filters by hard values
    |   |   |
    |   |   ├─ filter_applied  # FILTER column filled included those variants that didn't pass
    |   |   |
    |   |   └─ pass_only  # Only variants with FILTER=PASS
    |   |
    |   └─ VQSR  # (Option 2) Variant recalibration, filters by using truth and training data
    |       |
    |       ├─ model  # Models to be applied to VCFs
    |       |
    |       ├─ filter_applied  # FILTER column filled included those variants that didn't pass
    |       |
    |       └─ pass_only  # Only variants with FILTER=PASS
    |
    ├─* plink  # Contains plink files
    |
    ├─ genotypes
    |   |
    |   ├─ posteriors  # Calculate PP tag (posterior probabilites) and recaluclate MQ tag
    |   |
    |   ├─ filtered  # Set to missing genotypes with MQ < 20
    |   |
    |   └─ subsets  # Taking subsets of samples based on sequencing methods: GBS, WES, WGS
    |
    └─ haplotypes  # Phased/imputed VCFs
        |
        ├─ scaffolds  # Preliminary VCFs of genotypes incorportating pedigree information
        |
        └─ SHAPEIT4  # Phased/imputed VCFs using scaffold and original VCF


    The directories for relations from `relations.smk` are as follows. `plink` is created first. The others do not follow an order.
    However, these require the joint vcf file from `joint_call`:

    ├─ plink  # Contains plink files
    |
    ├─ relatedness  # Pairwise estimates of relatedness between samples
    |
    ├─ admixture  # Finds rates of admixture among samples
    |   |
    |   ├─ supervised  # Uses additional samples of known origin
    |   |
    |   └─ unsupervised  # No known origins
    |
    └─ aims  # Ancestry informative markers

    Directories for determining kinship using LASAR and SEEKIN:

    └─ kinship  # Allele frequencies and pairwise relatedness

    ''' > {output}
    """

rule recommended_resources:
    """Recommended directory layout for resources."""
    output: config["resources"] + "README.md"
    shell: """echo '''
    Recommended layout for resources:

    ├─ barcodes  # Table of barcodes for samples (used for GBS and AMP reads)
    |
    ├─ reads  # .fastq.gz file for each sample (R1 and R2)
    |
    ├─ ref_fna  # Reference genome
    |
    ├─ ref_vcf  # Reference .vcf.gz files
    |
    ├─ samples  # List of samples

    ''' > {output}
    """
