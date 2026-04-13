

# Find IBDs


# PRIMUS
# https://primus.gs.washington.edu/primusweb/res/documentation.html
# rule PRIMUS:
#     shell: """
#         run_PRIMUS.pl \
#         """


rule pixy:
    """Calculations using pixy. Want all sites included, even invariant ones to be most accurate.
    Especially need consideration when using something is not WGS. As we don't want to assume that those
    missing sites are invariant.
    """
    input:
        # Must be indexed vcf.gz
        vcfgz = config["results"] + "genotypes/pass/{dataset}.{subset}.SNP.autosomal.vcf.gz",
        tbi = config["results"] + "genotypes/pass/{dataset}.{subset}.SNP.autosomal.vcf.gz.tbi",
        populations = config["resources"] + "pop/MML_groups_from_Martha.fixed6.populations", # Headerless, tab-delimited, first column samples, second column population
        bed = "/master/abagwell/variant-analysis/results/rhesus/coverage/common_WES.cutoff=0.8.minDP=5.bed"
    output:
        config["results"] + "diversity/pixy/{dataset}.{subset}_pi.txt",
    params:
        out_dir = config["results"] + "diversity/pixy",
    threads: 10
    resources: nodes = 10
    conda: "../envs/admixture.yaml"
    shell: """
        pixy \
            --stats pi fst dxy \
            --vcf {input.vcfgz} \
            --populations <( \
                csvtk join {input.populations} <(bcftools query -l {input.vcfgz}) \
                    -t \
                    -f '1;1') \
            --window_size 10000 \
            --n_cores {threads} \
            --output_folder {params.out_dir} \
            --output_prefix {wildcards.dataset}.{wildcards.subset} \
            --bypass_invariant_check 'yes' \
            --bed_file {input.bed} \
        """


rule subset_BCF_to_population:
    """Split VCF into populations. Used for VCFtools nucleotide divergence."""
    input:
        bcf = config["results"] + "genotypes/pass/U42_WES.all2.SNP.autosomal.bcf",
        populations = config["resources"] + "pop/MML_groups_from_Martha.fixed6.populations",
    output:
        #bcf = config["results"] + "genotypes/subset/U42_WES.{}.SNP.autosomal.bcf",
        test = config["results"] + "genotypes/subset/split_populations.test.txt",
    conda: "../envs/common.yaml"
    shell: """
        while read POP; do \
            bcftools view {input.bcf} \
                -S '<( \
                    csvtk grep {input.populations} \
                        -t \
                        -f 2 \
                        -p $POP \
                    | csvtk cut \
                        -t \
                        -f 1 \
                    | paste -sd, \
                    )' \
            -Ob \
            -o /master/abagwell/variant-analysis/results/rhesus/genotypes/subset/U42_WES.$(POP).SNP.autosomal.bcf; \
        done <(csvtk cut -t -f 2 {input.populations} | sort | uniq | sed 's/ /_/g') \
        """


rule nucleotide_divergence_VCFtools:
    """Calculate nucleotide divergence with VCFtools."""
    input: 
        bcf = config["results"] + "genotypes/pass/{dataset}.{subset}.SNP.autosomal.bcf",
        csi = config["results"] + "genotypes/pass/{dataset}.{subset}.SNP.autosomal.bcf.csi",
    output:
     pi = config["results"] + "diversity/VCFtools/{dataset}.{subset}.sites.pi"
    conda: "../envs/common.yaml"
    shell: """
        vcftools \
            --bcf {input.bcf} \
            --site-pi \
            --stdout \
        > {output.pi} \
        """

rule inbreeding_coefficient:
    """Calculate inbreeding coefficient, which is a measure of heterozygosity per-individual."""
    input: vcf = config["results"] + "genotypes/pass/{dataset}.common_between_founding_cohorts.SNP.autosomal.bcf",
    output: het = config["results"] + "kinship/het/{dataset}.common_between_founding_cohorts.het",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools view {input.vcf} \
            --min-ac 1 \
            -Ov \
        | vcftools \
            --vcf - \
            --het \
            --stdout \
        > {output.het} \
        """
    # """
    #     bcftools view {input.vcf} \
    #         --min-ac 1 \
    #         -Ov
    #     | vcftools \
    #         --vcf - \
    #         --het \
    #         --stdout \
    #     | bcftools view \
    #         -Ob \
    #     > {output.het} \
    #     """

rule plink_het:
    """Computes, observed and expected autosomal homozygous genotype
    counts for each sample and method-of-moments."""
    input:
        # bed = config["results"] + "genotypes/pass/plink/{dataset}.{subset}.SNP.chr{chr}.bed",
        # bim = config["results"] + "genotypes/pass/plink/{dataset}.{subset}.SNP.chr{chr}.bim",
        # fam = config["results"] + "genotypes/pass/plink/{dataset}.{subset}.SNP.chr{chr}.fam",
        bcf = config["results"] + "genotypes/pass/{dataset}.{subset}.SNP.autosomal.bcf",
    output:
        het = config["results"] + "heterozygosity/PLINK/{dataset}.{subset}.SNP.autosomal.het"
    params:
        het = config["results"] + "heterozygosity/PLINK/{dataset}.{subset}.SNP.autosomal"
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        plink \
            --bcf {input.bcf} \
            --het \
            --out {params.het} \
            --allow-extra-chr \
        """
