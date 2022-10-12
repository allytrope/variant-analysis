"""Ancestry and pairwise relatedness estimation based on LASER and SEEKIN software, respectively.
LASER and SEEKIN-het are better than ADMIXTURE and lcMLkin for related individuals.
LASER is said to work on shotgun sequencing data, but requires a set of
reference individuals with "genome-wide SNP genotypes and ancestry information" available."""

CONFIG = config["kinship"]
DIR = config["results"] + "kinship/"
PREFIX = "no_missing_geno"

rule create_geno_and_site:
    """Create .geno and .site files for LASER."""
    input: ref_vcf = CONFIG["ref_vcf"],
           pop_ids = CONFIG["pop_ids"],
    output: geno = DIR + PREFIX + ".geno",
            site = DIR + PREFIX + ".site",
    params: prefix = DIR + PREFIX,
    shell: "vcf2geno \
                --inVcf {input.ref_vcf} \
                --out {params.prefix} \
                --updateID {input.pop_ids}"

rule create_bed:
    """Create .bed file, containing list of sites."""
    input: site = DIR + PREFIX + ".site",
    output: bed = DIR + PREFIX + ".bed",
    shell: "awk '{{if (NR>1) {{print $1,$2-1,$2,$3;}}}}' {input.site} > {output.bed}"

rule create_pileup:
    """Create .pileup file for a sample."""
    input: bed = DIR + PREFIX + ".bed",
           ref_fasta = config["ref_fasta"],
           ref_fasta_idx = config["ref_fasta"] + ".fai",
           bam = config["results"] + "alignments_recalibrated/{sample}.bam",  ## From bwa-mem in `variant_calling.smk`
    output: pileup = DIR + "pileup/{sample}.pileup",
    conda: "../envs/kinship.yaml"
    shell: "samtools sort {input.bam} | samtools mpileup \
                -q 30 \
                -Q 20 \
                -f {input.ref_fasta} \
                -l {input.bed} > {output.pileup}"

rule create_seq:
    """Create .seq file for LASER."""
    input: ref_fasta = config["ref_fasta"],
           site = DIR + PREFIX + ".site",
           pileups = expand("{dir}pileup/{sample}.pileup", dir=DIR, sample=SAMPLE_NAMES)
    output: seq = DIR + PREFIX + ".seq",
    # SET IN PATH OR MOVE
    # `-b` and `-i` flags are optional
    shell: "python pileup2seq.py {input.pileups} \
                -f {input.ref_fasta} \
                -m {input.site} \
                -o {output.seq}"

rule create_coord:
    """Create .coord file for LASER by running in PCA mode."""
    input: geno = DIR + PREFIX + ".geno",
    output: coord = DIR + PREFIX + ".RefPC.coord",
            var = DIR + PREFIX + ".RefPC.var",
            grm = DIR + PREFIX + ".RefPC.grm",  # Only for when `-pca 1`. When `3`, makes `.RefPC.load` instead.
    params: pca_mode = 1,
            prefix = DIR + PREFIX,
            conf = DIR + "laser.conf",
    shell: "laser \
                -pca {params.pca_mode} \
                -o {params.prefix} \
                -p {params.conf}"

rule laser:
    """Estimate individual ancestry background."""
    input: geno = DIR + PREFIX + ".geno",
           seq = DIR + PREFIX + ".seq",
           ref_coord = DIR + PREFIX + ".RefPC.coord",
    output: seq_coord = DIR + PREFIX + ".SeqPC.coord",
            #ind_cov = DIR + OUT_PREFIX + ".ind.cov",  # When `-cov 1`
            #loc_cov = DIR + OUT_PREFIX + ".loc.cov",  # When `-cov 1`
    params: prefix = DIR + PREFIX,
            conf = DIR + "laser.conf",  # parameterfile
    shell: "laser \
                -g {input.geno} \
                -s {input.seq} \
                -c {input.ref_coord} \
                -o {params.prefix} \
                -p {params.conf}"

rule model_AF:
    """Model allele frequencies as linear functions of PCs based on the ancestry reference panel."""
    input: vcf = config["ref_vcf"],  # genotypes of reference individuals (.vcf.gz)
           ref_coord = DIR + "RefPC.coord",  # PCA coordinates of reference individuals
    output: DIR + "AF.model",
    threads: 1
    params: k = CONFIG["dim"],  # Number of PCs used to model AF. Default is 2.
    shell: "seekin modelAF \
                -i {input.vcf} \
                -c {input.ref_coord} \
                -k {params.k} \
                -o {output}"

rule get_AF:
    """Estimate individual-specific allele frequencies of study individuals."""
    input: ref_coord = DIR + PREFIX + ".RefPC.coord",  # coordinates of study individuals in the reference space
           AFmodel = DIR + "AF.model",
    output: indivAF = DIR + "indivAF.vcf.gz",
    threads: 1
    params: k = CONFIG["dim"], # Number of PCs used to compute allele frequencies.
    shell: "seekin getAF \
                -i {input.ref_coord} \
                -b {input.AFmodel} \
                -k {params.k} \
                -o {output.indivAF}"

rule pairwise_kinship:
    """Estimate relatedness with SEEKIN-het for related individuals and admixture."""
    input: indivAF = DIR + PREFIX + ".indivAF.vcf.gz",
           ref_vcf = CONFIG["ref_vcf"],  ## VCF of genotypes or dosages of study individuals
    output: log = DIR + PREFIX + ".het.log",
            kin = DIR + PREFIX + ".het.kin",
            inbreed = DIR + PREFIX + ".het.inbreed",
            matrix = DIR + PREFIX + ".het.matrix",
            matrixID = DIR + PREFIX + ".het.matrixID",
    threads: 10  # Default 10
    params: prefix = DIR + PREFIX + ".het",
    shell: "seekin kinship \
            -i {input.ref_vcf} \
            -f {input.indivAF} \
            -p het \
            -t {threads} \
            -o {params.prefix}"

