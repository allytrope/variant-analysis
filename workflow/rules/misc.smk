"""Micellaneous rules."""

rule concat_chromosomes:
    """Concatenate chromosomes."""
    input:
        vcfs = lambda wildcards: expand(config["results"] + "{path}/{dataset}.{mode}.chr{chr}.vcf.gz",
            path=wildcards.path,
            dataset=wildcards.dataset,
            mode=wildcards.mode,
            chr=CHROMOSOMES),
        chromosomes = config["results"] + "db/chromosomes.list",
    output:
        vcf = config["results"] + "{path}/{dataset}.{mode}.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools concat {input.vcfs} \
            -o {output.vcf} \
            -Oz \
        """