
rule kin_akt:
    """Calculate kinship using the tool akt."""
    input:
        # vcf = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.{mode}.vcf.gz",
        # tbi = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.{mode}.vcf.gz.tbi",
        # regions = config["results"] + "relatedness/pca/SHAPEIT5_WGS/{dataset}.{mode}.pca_regions.bed",
        #bcf = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.bcf",
        bcfs = expand(config["results"] + "genotypes/pass/{{dataset}}.{{subset}}.{{mode}}.chr{chr}.bcf",
            chr=AUTOSOMES),
    output:
        txt = config["results"] + "kinship/akt/{dataset}.{subset}.{mode}.kin.txt",
    threads: 24
    resources: nodes = 24
    conda: "../envs/kinship.yaml"
    shell: """
        bcftools concat {input.bcfs} -Ob \
        | akt kin \
            --force \
            --threads {threads} \
        > {output.txt} \
        """

# rule merge_plink:
#     input:
#         beds = lambda wildcards: expand(config["results"] + "genotypes/plink/{dataset}.{subset}.SNP.chr{chr}.bed",
#             dataset=wildcards.dataset,
#             subset=wildcards.subset,
#             chr=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"]),
#     output:
#         tmp = config["results"] + "genotypes/plink/{dataset}.{subset}.SNP.autosomal.list",
#         merged = config["results"] + "genotypes/plink/{dataset}.{subset}.SNP.autosomal.bed",
#     params: 
#         first_bed = config["results"] + "genotypes/plink/{dataset}.{subset}.SNP.chr1",
#         other_beds = lambda wildcards: expand(config["results"] + "genotypes/plink/{dataset}.{subset}.SNP.chr{chr}",
#             dataset=wildcards.dataset,
#             subset=wildcards.subset,
#             chr=["2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"]),
#         merged = config["results"] + "genotypes/plink/{dataset}.{subset}.SNP.autosomal",
#     conda: "../envs/rvtests.yaml"
#     shell: """
#         for i in {params.other_beds}; do echo $i >> {output.tmp}; done; \
#         plink --file {params.first_bed} --merge-list {output.tmp} --make-bed --out {params.merged} \
#         """

rule estimate_relatedness_king:
    """Estimate relatedness between individuals. It is recommended by KING not to prune or filter any "good" SNPs for this."""
    input: 
        bed = config["results"] + "genotypes/pass/plink/{dataset}.{subset}.SNP.autosomal.bed",
        # bim = config["results"] + "genotypes/plink/{dataset}.{subset}.SNP.autosomal.bim",
        # fam = config["results"] + "genotypes/plink/{dataset}.{subset}.SNP.autosomal.fam",
    output:
        kinship = config["results"] + "kinship/KING/{dataset}.{subset}.SNP.autosomal.kin",
    params:
        prefix = config["results"] + "kinship/KING/{dataset}.{subset}.SNP.autosomal",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        king \
            -b {input.bed} \
            --kinship \
            --prefix {params.prefix} \
        """




rule estimate_relatedness_king_per_chr:
    """Estimate relatedness between individuals."""
    input: 
        # bed = config["results"] + "kinship/plink/GBS_WES_WGS.chr{chr}.same_fid.bed",
        # bim = config["results"] + "kinship/plink/GBS_WES_WGS.chr{chr}.same_fid.bim",
        # fam = config["results"] + "kinship/plink/GBS_WES_WGS.chr{chr}.same_fid.fam",
        bed = config["results"] + "genotypes/pass/plink/{dataset}.{subset}.SNP.chr{chr}.bed",
        bim = config["results"] + "genotypes/pass/plink/{dataset}.{subset}.SNP.chr{chr}.bim",
        fam = config["results"] + "genotypes/pass/plink/{dataset}.{subset}.SNP.chr{chr}.fam",
    output:
        kinship = config["results"] + "kinship/KING/{dataset}.{subset}.SNP.chr{chr}.kin",
    params:
        prefix = config["results"] + "kinship/KING/{dataset}.{subset}.SNP.chr{chr}",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        king \
            -b {input.bed} \
            --kinship \
            --prefix {params.prefix} \
        """

rule king_matrix_inter:
    """Intermediate step."""
    input:
        kinship = config["results"] + "kinship/KING/{dataset}.{subset}.SNP.autosomal.kin",
    output:
        matrix = temp(config["results"] + "kinship/KING/{dataset}.{subset}.SNP.autosomal.inter.matrix"),
    run:
        import pandas as pd
        df = pd.read_csv(input.kinship, delimiter='\t')

        # make unique, sorted, common index
        idx = sorted(set(df['ID1']).union(df['ID2']))

        # reshape
        (df.pivot(index='ID1', columns='ID2', values='Kinship')
        .reindex(index=idx, columns=idx)
        .fillna(0, downcast='infer')  # Change NAs to 0
        .pipe(lambda x: x+x.values.T)  # Fill in transposed side of table
        ).to_csv(output.matrix, sep='\t')

# rule king_matrix:
#     """Set diagonals to 0.5."""
#     input:
#         matrix = config["results"] + "kinship/KING/{dataset}.{subset}.SNP.autosomal.inter.matrix",
#     output:
#         matrix = config["results"] + "kinship/KING/{dataset}.{subset}.SNP.autosomal.matrix",
#     shell: """
#         awk '{{$NR=0.5; print}}' input.matrix \
#         > {output.matrix} \
#         """

rule king_matrix_PMx:
    """Convert KING output into a matrix format required by PMx."""
    input:
        matrix = config["results"] + "kinship/KING/{dataset}.{subset}.SNP.autosomal.inter.matrix",
    output:
        PMx_matrix = config["results"] + "kinship/KING/{dataset}.{subset}.SNP.autosomal.PMx.matrix",
    shell: """
        cut {input.matrix} -f 1 | cut -c 4- | cut -d _ -f 1 \
        | sed '1d' \
        >> {output.PMx_matrix}; \
        cut {input.matrix} \
            -f 2- \
        | sed '1d' \
        | awk '{{$NR=0.5; print}}' \
        >> {output.PMx_matrix}; \
        """
        

