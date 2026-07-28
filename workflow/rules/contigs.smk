## List chromosomes or subsets

rule list_contigs:
    """Find all contigs from headers of reference genome"""
    input:
        ref = config["ref_fasta"],
    output:
        config["resources"] + "ref_fna/contigs.list",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        zcat {input.ref} | grep "^>" | cut -c 2- | cut -d " " -f 1 > {output}
        """

rule list_chromosomes:
    """Find all chromosomes.
    
     Each line of the output file is just one chromosome's name.
     The search only keeps numbered chromosomes (not those prefixed with "chr" or any other letters)
     as well as X, Y, and MT. Unplaced contigs are ignored."""
    input:
        contigs = config["resources"] + "ref_fna/contigs.list",
    output:
        config["resources"] + "ref_fna/chromosomes.list",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        cat {input.contigs} | grep -x -E "^[0-9]+|X|Y|MT" > {output}
        """

rule list_autosomes:
    """Find all autosomes. That is, only the numbered chromosomes. Does not include unplaced contigs."""
    input:
        chromosomes = config["resources"] + "ref_fna/chromosomes.list",
    output:
        config["resources"] + "ref_fna/autosomes.list",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        cat {input.chromosomes} | grep -vxE "X|Y|MT" > {output}
        """

## BAMs
# rule split_bam_chromosomes:


## VCFs

rule concat_chromosomes:
    """Concatenate chromosomes."""
    input:
        vcfs = expand(config["results"] + "{{path}}/{{dataset}}.{{mode}}.chr{chr}.vcf.gz",
            chr=CHROMOSOMES),
        chromosomes = config["resources"] + "ref_fna/chromosomes.list",
    output:
        vcf = config["results"] + "{path}/{dataset}.{mode}.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools concat {input.vcfs} \
            -o {output.vcf} \
            -Oz \
        """

rule bcfs_to_autosomal_bcf:
    """Concatenate autosomes."""
    wildcard_constraints:
        ext = "bcf",
    input:
        bcfs = expand(config["results"] + "{{path}}.chr{chr}.{{ext}}",
            chr=AUTOSOMES),
        chromosomes = config["resources"] + "ref_fna/chromosomes.list",
    output:
        bcf = config["results"] + "{path}.autosomal.{ext}",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools concat {input.bcfs} \
            -o {output.bcf} \
            -Ob \
        """

# rule vcfgzs_to_autosomal_vcfgz:
#     """Concatenate autosomes."""
#     wildcard_constraints:
#         ext = "vcf.gz",
#     input:
#         vcfs = lambda wildcards: expand(config["results"] + "{path}.chr{chr}.{ext}",
#             path=wildcards.path,
#             chr=[chrom for chrom in CHROMOSOMES if chrom not in ["X", "Y", "MT"]],
#             ext=wildcards.ext),
#         chromosomes = config["resources"] + "ref_fna/chromosomes.list",
#     output:
#         vcf = config["results"] + "{path}.autosomal.{ext}",
#     threads: 1
#     resources: nodes = 1
#     conda: "../envs/common.yaml"
#     shell: """
#         bcftools concat {input.vcfs} \
#             -o {output.vcf} \
#             -Oz \
#         """

rule bcfs_to_autosomal_vcfgz:
    """Concatenate autosomes."""
    input:
        bcfs = lambda wildcards: expand(config["results"] + "{path}.chr{chr}.bcf",
            path=wildcards.path,
            chr=[chrom for chrom in CHROMOSOMES if chrom not in ["X", "Y", "MT"]],
            #ext=wildcards.ext
        ),
        chromosomes = config["resources"] + "ref_fna/chromosomes.list",
    output:
        vcf = config["results"] + "{path}.autosomal.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools concat {input.bcfs} \
            -o {output.vcf} \
            -Oz \
        """

#ruleorder: bcfs_to_autosomal_vcfgz > vcfgzs_to_autosomal_vcfgz

# Ambiguous with biallelics_by_mode
rule split_chromosomes:
    """Split into chromosomes."""
    input:
        vcf = config["results"] + "{path}/{name}.bcf",
        csi = config["results"] + "{path}/{name}.bcf.csi",
    output:
        vcf = config["results"] + "{path}/{name}.chr{chr}.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools view {input.vcf} \
            -r {wildcards.chr} \
            -Ob \
            -o {output.vcf}
        """
