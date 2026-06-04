"""Rules related to principle component analysis."""


rule pca_regions:
    """Select regions for PCA at random."""
    input:
        vcf = config["results"] + "haplotypes/SHAPEIT4/with_scaffold/{dataset}.{mode}.vcf.gz",
    output:
        regions = config["results"] + "relatedness/pca/{dataset}.{mode}.pca_regions.bed",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools view {input.vcf} -H \
        | shuf -n 8000 \
        | cut -f 1,2 \
        | sort -n -k 2 \
        | awk 'BEGIN {{OFS="\t"}} {{print $1,$2-1,$2}}' \
        > {output.regions} \
        """

rule pca_regions_VCF:
    """Select regions for PCA at random as VCF."""
    input:
        vcf = config["results"] + "haplotypes/SHAPEIT4/with_scaffold/{dataset}.{mode}.vcf.gz",
    output:
        regions = config["results"] + "relatedness/pca/{dataset}.{mode}.pca_regions.sorted.vcf.gz",
        tmp = temp(config["results"] + "relatedness/pca/{dataset}.{mode}.pca_regions.unsorted.vcf"),
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools view {input.vcf} -hG > {output.tmp}; \

        bcftools annotate -x INFO,FORMAT {input.vcf} \
        | bcftools view -GH \
        | shuf -n 8000 \
        >> {output.tmp}; \

        bcftools sort {output.tmp} \
        | bcftools view -Oz \
        > {output.regions} \
        """

rule pca_akt:
    """PCA using the tool akt."""
    input:
        vcf = config["results"] + "haplotypes/SHAPEIT5/{dataset}.{mode}.vcf.gz",
        tbi = config["results"] + "haplotypes/SHAPEIT5/{dataset}.{mode}.vcf.gz.tbi",
        #regions = config["results"] + "relatedness/pca/{dataset}.{mode}.pca_regions.bed.gz",
        regions = config["results"] + "relatedness/pca/{dataset}.{mode}.pca_regions.sorted.vcf.gz",
        regions_tbi = config["results"] + "relatedness/pca/{dataset}.{mode}.pca_regions.sorted.vcf.gz.tbi",
    output:
        vcf = config["results"] + "relatedness/pca/{dataset}.{mode}.pca.{npca}.vcf.gz",
        txt = temp(config["results"] + "relatedness/pca/{dataset}.{mode}.pca.{npca}.txt"),
    threads: 1
    resources: nodes = 1
    conda: "../envs/kinship.yaml"
    shell: """
        akt pca {input.vcf} \
            -R {input.regions} \
            --npca {wildcards.npca} \
            -Oz \
            -o {output.vcf} \
        > {output.txt} \
        """

rule pca_to_tsv:
    """Reformat PCA file into TSV."""
    input:
        txt = config["results"] + "relatedness/pca/{dataset}.{mode}.pca.{npca}.txt",
    output:
        tsv = config["results"] + "relatedness/pca/{dataset}.{mode}.pca.{npca}.tsv",
    # There is likely a way to use fewer substitutions, but this works as is.
    shell:
        """
        sed 's/\\t/ /g;s/   /\\t/g;s/  / /g;s/ /\\t/g;s/\\t\\t/\\t/g' {input.txt} > {output.tsv}
        """

## PLINK PCA ##
rule plink_linkage_pruning:
    shell: """
        plink \
            --vcf $VCF \
            --double-id \
            --allow-extra-chr \
            --set-missing-var-ids @:# \
            --indep-pairwise 50 10 0.1 --out cichlids \
        """

rule plink_pca:
    shell: """
        plink \
            --vcf $VCF \
            --double-id \
            --allow-extra-chr \
            --set-missing-var-ids @:# \
            --extract cichlids.prune.in \
            --make-bed \
            --pca \
            --out cichlids \
        """