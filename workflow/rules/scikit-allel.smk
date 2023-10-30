"""Rules for calculating stats using scikit-allel."""

CHROMOSOMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X"]

## General rules used by others

rule subpop_vcf:
    """Create a VCF containing only individuals of one subpopulation. Only works if all samples are actually in VCF."""
    # wildcard_constraints:
    #     chrom = "^chr|\."
    input:
        vcf = config["results"] + "{path}/{description}.chr{chr}.bcf",
        samples = lambda wildcards: config["pops"][int(wildcards.int)],
    output:
        vcf = config["results"] + "{path}/{description}.pop{int}.chr{chr}.vcf.gz",
    conda: "../envs/bio.yaml"
    shell: """
        bcftools view {input.vcf} \
            -S {input.samples} \
            -Oz \
            -o {output.vcf} \
        """

# This rule doesn't work with the changes, but not sure why.
rule concat_chrom_chromosomes:
    """Merge chromosomal stats from scikit_allel jobs."""
    wildcard_constraints:
        #statistic = "divergence|diversity|expected_heterozygostiy",
        #seq = "WGS",
    input:
        # pickles = lambda wildcards: expand(config["results"] + "scikit-allel/{statistic}/{seq}/{dataset}.{mode}.{window_size}kbp.chr{chr}.pickle", 
        #     statistic=wildcards.statistic,
        #     seq=wildcards.seq,
        #     dataset=wildcards.dataset,
        #     mode=wildcards.mode,
        #     window_size=wildcards.window_size,
        #     chr=CHROMOSOMES),
        pickles = lambda wildcards: expand(config["results"] + "scikit-allel/{path}/{description}.chr{chr}.pickle", 
            path=wildcards.path,
            description=wildcards.description,
            chr=CHROMOSOMES),
    output:
        #pickle = config["results"] + "scikit-allel/{statistic}/{seq}/{dataset}.{mode}.{window_size}kbp.pickle"
        pickle = config["results"] + "scikit-allel/{path}/{description}.all.pickle",
    conda: "../envs/scikit.yaml"
    script:
        "../scripts/scikit-allel/concat_chrom_stats.py"

rule concat_pop_stats:
    """Merge subpopulations stats from scikit_allel jobs."""
    input:
        pickles = lambda wildcards: expand(config["results"] + "scikit-allel/{path}/{description}.pop{subpop}.all.pickle", 
            path=wildcards.path,
            description=wildcards.description,
            subpop=range(0, len(config["pops"]))),
    output:
        pickle = config["results"] + "scikit-allel/{path}/{description}.merged.pickle",
    conda: "../envs/scikit.yaml"
    script:
        "../scripts/scikit-allel/concat_pop_stats.py"


# rule concat_allel_chromosomes:
#     """Merge chromosomal stats from scikit_allel jobs."""
#     wildcard_constraints:
#         #statistic = "divergence|diversity|expected_heterozygostiy",
#         #seq = "WGS",
#     input:
#         pickles = lambda wildcards: expand(config["results"] + "scikit-allel/expected_heterozygosity/WGS/SNPRC_WGS_WES.SNP.chr{chr}.pickle", 
#             statistic=wildcards.statistic,
#             seq=wildcards.seq,
#             dataset=wildcards.dataset,
#             mode=wildcards.mode,
#             chr=CHROMOSOMES),
#         # pickles = lambda wildcards: expand(config["results"] + "scikit-allel/{path}/{description}.chr{chr}.pickle", 
#         #     path=wildcards.path,
#         #     description=wildcards.description,
#         #     chr=CHROMOSOMES),
#     output:
#         pickle = config["results"] + "scikit-allel/expected_heterozygosity/WGS/SNPRC_WGS_WES.SNP.pickle"
#     conda: "../envs/scikit.yaml"
#     script:
#         "../scripts/scikit-allel/concat_chrom_stats.py"

rule concat_allel_chromosomes_by_window:
    """Merge windowed stats from scikit_allel jobs. Used after merging by chromosome."""
    input:
        pickles = lambda wildcards: expand(config["results"] + "scikit-allel/{statistic}/{seq}/{dataset}.{mode}.{window_size}kbp.pickle", 
            statistic=wildcards.statistic,
            seq=wildcards.seq,
            dataset=wildcards.dataset,
            mode=wildcards.mode,
            window_size=config["window_sizes"]),
    output:
        pickle = config["results"] + "scikit-allel/{statistic}/{seq}/{dataset}.{mode}.pickle"
    conda: "../envs/scikit.yaml"
    script:
        "../scripts/scikit-allel/concat_chrom_stats_by_window.py"

## Divergence

rule divergence:
    input:
        vcf = config["results"] + "haplotypes/SHAPEIT5_{seq}/{dataset}.{mode}.chr{chr}.vcf.gz",
        pops = config["pops"],
    output:
        pickle = temp(config["results"] + "scikit-allel/divergence/{seq}/{dataset}.{mode}.{window_size}kbp.chr{chr}.pickle"),
    conda: "../envs/scikit.yaml"
    script:
        "../scripts/scikit-allel/divergence.py"

rule render_divergence:
    """Create Altair plot in html."""
    input:
        pickle = lambda wildcards: expand(config["results"] + "scikit-allel/divergence/{seq}/{dataset}.{mode}.{window_size}kbp.pickle",
            seq=wildcards.seq,
            dataset=wildcards.dataset,
            mode=wildcards.mode,
            window_size=wildcards.window_size)
    output:
        html = config["results"] + "scikit-allel/divergence/{seq}/{dataset}.{mode}.{window_size}kbp.html",
    conda: "../envs/graph.yaml"
    # notebook:
    #     "../notebooks/scikit-allel/render_divergence.py.ipynb"
    script:
        "../notebooks/scikit-allel/render_divergence.py"

rule render_divergence_all_windows:
    """Create Altair plot in html."""
    input:
        pickles = lambda wildcards: expand(config["results"] + "scikit-allel/divergence/{seq}/{dataset}.{mode}.pickle",
            seq=wildcards.seq,
            dataset=wildcards.dataset,
            mode=wildcards.mode,)
    output:
        html = config["results"] + "scikit-allel/divergence/{seq}/{dataset}.{mode}.html",
    conda: "../envs/graph.yaml"
    notebook:
        "../notebooks/scikit-allel/render_divergence_all_windows.py.ipynb"

## Diversity

rule diversity:
    input:
        config["results"] + "haplotypes/SHAPEIT5_WGS/subpop/{subpop}/{dataset}.{mode}.chr{chr}.vcf.gz",
    output:
    conda: "../envs/scikit.yaml"
    script:
        "../scripts/scikit-allel/diversity.py"

## Expected heterozygosity

rule expected_heterozygosity:
    """Expected rate of heterozygosity per variant under Hardy-Weinberg equilibrium."""
    input:
        #vcf = config["results"] + "haplotypes/SHAPEIT5_{seq}/{dataset}.{mode}.pop{subpop}.chr{chr}.vcf.gz",
        vcf = config["results"] + "haplotypes/SHAPEIT5_{seq}/{dataset}.{mode}.chr{chr}.vcf.gz",
        subpop = lambda wildcards: config["pops"][int(wildcards.subpop)],
    output:
        pickle = temp(config["results"] + "scikit-allel/expected_heterozygosity/{seq}/{dataset}.{mode}.pop{subpop}.chr{chr}.pickle"),
    conda: "../envs/scikit.yaml"
    script:
        "../scripts/scikit-allel/expected_heterozygosity.py"

# rule render_expected_heterozygosity:
#     input:
#         pickle = config["results"] + "scikit-allel/expected_heterozygosity/{seq}/{dataset}.{mode}.chr{chr}.pickle",
#     output:
#     conda: "../envs/graph.yaml"
#     script:
#         "../notebooks/scikit-allel/expected_heterozygosity.py.ipynb"

# rule sgkit_TajimaD:
#     """Es"""
#     input:
#     output:
#     conda: "../envs/scikit.yaml"
#     run: """
        
#         """

# rule scikit_allel_delta_TajimaD:
#     """Difference in Tajima's D between two populations with moving windows."""
#     input:
#     output:
#     conda: "../envs/scikit.yaml"
#     run: """
#         import allel
#         allel.moving_delta_tajima_d()
#         """

# rule scikit_allel_XPEHH:
#     """Unstandardized cross-population extended haplotype homozygosity."""
#     input:
#     output:
#     conda: "../envs/scikit.yaml"
#     shell: """
#         import allel
#         allel.xpehh()
#         """



