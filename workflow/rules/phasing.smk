"""Contain rules for genotype refinement and haplotype estimation.
This comes after hard_filter.smk or modified to come after variant_recalibration.smk."""

CONFIG = config["phasing"]


## Genotype Refinement

rule genotype_posteriors:
    """Calculate genotype posterior probabilties. This adds a PP field in the FORMAT column of the VCF for the genotype posteriors.
    Also, the GQ and GT fields may also be subject to change. These are caluclated using information from the other samples in the file
    and also specifically trio information."""
    input:
        vcf = config["results"] + "hard_filtered/pass_only/biallelic_split/{workspace}.{mode}.chr{chr}.vcf.gz",
        vcf_idx = config["results"] + "hard_filtered/pass_only/biallelic_split/{workspace}.{mode}.chr{chr}.vcf.gz.tbi",
        ped = config["resources"] + "pedigree/trios.ped",
    output:
        config["results"] + "genotypes/posteriors/{workspace}.{mode}.chr{chr}.posteriors.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    shell: """
        gatk --java-options "-Xmx8g" CalculateGenotypePosteriors \
        -V {input.vcf} \
        -ped {input.ped} \
        -O {output} \
        --tmp-dir ~/tmp/{rule}"""

# Already have hard_filter.smk rules
rule genotype_filtration:
    """Filtering genotypes for GQ < 20 (which represents a 99% chance of being correct).
    This adds the filter tag if fails and sets "./." as new genotype."""
    input:
        vcf = config["results"] + "genotypes/posteriors/{workspace}.{mode}.chr{chr}.posteriors.vcf.gz",
        vcf_idx = config["results"] + "genotypes/posteriors/{workspace}.{mode}.chr{chr}.posteriors.vcf.gz.tbi",
    output:
        config["results"] + "genotypes/filtered/{workspace}.{mode}.chr{chr}.filtered.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    shell: """
        gatk --java-options "-Xmx8g" VariantFiltration \
            -V {input.vcf} \
            --genotype-filter-expression "GQ < 20" \
            --genotype-filter-name "GQ20" \
            --set-filtered-genotype-to-no-call \
            -O {output}"""

# rule mark_de_novo:
#     """Tag possible de novo mutations."""
#     input:
#     output:
#     shell: """
#         gatk"""

rule split_by_primary_sequencing_type:
    """Split into GBS, WES, WGS."""
    wildcard_constraints: subset = "GBS|WES|WGS"
    input:
        vcf = config["results"] + "genotypes/filtered/GBS_WES_WGS.{mode}.chr{chr}.filtered.vcf.gz",
        samples = config["results"] + "genotypes/{subset}.list",  # Must be manually made
    output:
        config["results"] + "genotypes/subsets/{subset}/{subset}.{mode}.chr{chr}.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools view {input.vcf} \
            -S {input.samples} \
            --min-ac 1 \
            -Oz \
            -o {output}"""

rule split_by_multiple_sequencing_types:
    """Create GBS_WGS or WES_WGS."""
    wildcard_constraints: subset = "GBS|WES"
    input:
        vcf = config["results"] + "genotypes/filtered/GBS_WES_WGS.{mode}.chr{chr}.filtered.vcf.gz",
        samples = config["results"] + "genotypes/{subset}_WGS.list",  # Must be manually made
        targets = config["results"] + "genotypes/subsets/{subset}/{subset}.{mode}.chr{chr}.vcf.gz",
    output:
        config["results"] + "genotypes/subsets/{subset}_WGS/{subset}_WGS.{mode}.chr{chr}.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools view {input.vcf} \
            --samples-file {input.samples} \
            --targets-file {input.targets} \
            --min-ac 1 \
            -Oz \
            -o {output}"""


### Phasing

# Not currently being used
rule whatshap:
    """Haplotype assembly to be used as phase set for SHAPEIT4."""
    input: 
        vcf = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.filtered.vcf.gz",  # Taking hard filtered .vcf
        bams = expand(config["results"] + "alignments_recalibrated/{sample}.bam", sample=SAMPLE_NAMES),  # Should include all samples in VCF
        ref_fasta = config["ref_fasta"],
        ped = ""  # PLINK PED file (but only uses columns 2, 3, and 4)
    output:
        phased = config["results"] + "haplotypes/whatshap/{dataset}/{dataset}.{mode}.chr{chr}.phased.vcf.gz",
    threads: 1
    resources: nodes = 1
    shell: """
        whatshap phase {input.vcf} {input.bams} \
            --reference={input.ref_fasta} \
            --ped {input.ped} \
            -o {output.vcf} \
        """

rule make_scaffold:
    input:
        vcf = config["results"] + "genotypes/subsets/{dataset}/{dataset}.SNP.chr{chr}.vcf.gz",  # Taking hard filtered .vcf
        vcf_idx = config["results"] + "genotypes/subsets/{dataset}/{dataset}.SNP.chr{chr}.vcf.gz.tbi",
        fam = config["resources"] + "pedigree/parents.txt",
    output:
        scaffold = config["results"] + "haplotypes/scaffolds/{dataset}/{dataset}.{mode}.chr{chr}.scaffold.vcf.gz",
    threads: 1
    resources: nodes = 1
    shell: """
        makeScaffold \
            --gen {input.vcf} \
            --fam {input.fam} \
            --reg {wildcards.chr} \
            --out {output.scaffold} \
        """

rule shapeit4:
    input:
        vcf = config["results"] + "genotypes/subsets/{dataset}/{dataset}.{mode}.chr{chr}.vcf.gz",  # Taking hard filtered .vcf
        vcf_idx = config["results"] + "genotypes/subsets/{dataset}/{dataset}.{mode}.chr{chr}.vcf.gz.csi",
        scaffold = config["results"] + "haplotypes/scaffolds/{dataset}/{dataset}.{mode}.chr{chr}.scaffold.vcf.gz",
        scaffold_idx = config["results"] + "haplotypes/scaffolds/{dataset}/{dataset}.{mode}.chr{chr}.scaffold.vcf.gz.csi",
        reference = config["results"] + "haplotypes/scaffolds/WGS/WGS.SNP.chr{chr}.scaffold.vcf.gz",
        reference_idx = config["results"] + "haplotypes/scaffolds/WGS/WGS.SNP.chr{chr}.scaffold.vcf.gz.csi",
    output:
        phased = config["results"] + "haplotypes/SHAPEIT4/{dataset}/{dataset}.{mode}.chr{chr}.phased.vcf.gz",
    # When using SHAPIT4.1 or greater, --map not required (though surely still helpful).
    # --chr appears to still be mandatory even if there is only one chromosome in file.
    # --scffold is useful for incorporation pedigree information.
    log: config["results"] + "haplotypes/SHAPEIT4/log/{dataset}.{mode}.chr{chr}.log",
    threads: 12
    resources: nodes = 12
    conda: "../envs/shapeit4.yaml"
    shell: """
        shapeit4 \
            --input {input.vcf} \
            --reference {input.reference} \
            --scaffold {input.scaffold} \
            --region {wildcards.chr} \
            --sequencing \
            --output {output.phased} \
            --log {log} \
            --thread {threads} \
        """

