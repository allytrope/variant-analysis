# TODO: Generally apply to n-datasets
rule left_join_VCFs:
    """Merge VCF datasets, keeping additional samples for an individual from left dataset only if in both."""
    input:
        vcf1 = config["results"] + "genotypes/pass/{dataset1}.all.{mode}.chr{chr}.bcf",
        csi1 = config["results"] + "genotypes/pass/{dataset1}.all.{mode}.chr{chr}.bcf.csi",
        vcf2 = config["results"] + "genotypes/pass/{dataset2}.all.{mode}.chr{chr}.bcf",
        csi2 = config["results"] + "genotypes/pass/{dataset2}.all.{mode}.chr{chr}.bcf.csi",
    params:
        regions = config["results"] + "genotypes/pass/{dataset1}.all.{mode}.chr{chr}.regions.tsv",
    output:
        temp = temp(config["results"] + "genotypes/pass/{dataset1}+{dataset2}.temp.{mode}.chr{chr}.bcf"),
        vcf = config["results"] + "genotypes/pass/{dataset1}+{dataset2}_left_join.all.{mode}.chr{chr}.bcf",
    shell: """
        bcftools view {input.vcf1} -H | cut -f 1,2 > {params.regions};
        bcftools merge {input.vcf1} {input.vcf2} \
            --force-samples \
            -R {params.regions} \
            -Ob \
        > {output.temp};
        bcftools view {output.temp} \
            -S <(bcftools query -l {output.temp} | grep -v -F '^2:') \
        > {output.vcf};
    """

rule outer_join_samples:
    """Merge VCF datasets, keeping all samples from both datasets."""
    input:
        vcf1 = config["results"] + "genotypes/filtered/{dataset1}.{mode}.chr{chr}.vcf.gz",
        csi1 = config["results"] + "genotypes/filtered/{dataset1}.{mode}.chr{chr}.vcf.gz.tbi",
        vcf2 = config["results"] + "genotypes/filtered/{dataset2}.{mode}.chr{chr}.vcf.gz",
        csi2 = config["results"] + "genotypes/filtered/{dataset2}.{mode}.chr{chr}.vcf.gz.tbi",
    output:
        vcf = config["results"] + "genotypes/filtered/{dataset1}+{dataset2}_outer_join.{mode}.chr{chr}.vcf.gz",
    shell: """
        bcftools merge {input.vcf1} {input.vcf2} \
            --force-samples \
            -Oz \
        > {output.vcf};
    """