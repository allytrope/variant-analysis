"""Contain rules for hard filtering. Runs after variant_calling.smk and as an alternative to variant_recalibration.smk."""

CONFIG = config["hard_filter"]


# def filters(mode):
#     """Use in rule hard_filter."""
#     if mode == "SNP" or mode == "both":
#         return f"""-filter "QUAL < {CONFIG["QUAL"]}" --filter-name "QUAL{CONFIG["QUAL"]}" \
#                   -filter "QD < {CONFIG["QD"]}" --filter-name "QD{CONFIG["QD"]}" \
#                   -filter "SOR > {CONFIG["SOR"]}" --filter-name "SOR{CONFIG["SOR"]}" \
#                   -filter "FS > {CONFIG["FS"]}" --filter-name "FS{CONFIG["FS"]}" \
#                   -filter "MQ < {CONFIG["MQ"]}" --filter-name "MQ{CONFIG["MQ"]}" \
#                   -filter "MQRankSum < {CONFIG["MQRankSum"]}" --filter-name "MQRankSum{CONFIG["MQRankSum"]}" \
#                   -filter "ReadPosRankSum < {CONFIG["ReadPosRankSum"]}" --filter-name "ReadPosRankSum{CONFIG["ReadPosRankSum"]}" """
#     elif mode == "indel":
#         return """-filter "QUAL < 30.0" --filter-name "QUAL30" \
#                   -filter "QD < 2.0" --filter-name "QD2" \
#                   -filter "FS > 200.0" --filter-name "FS200" \
#                   -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" """
rule pass_only_hard_filter:
    """Remove variants by filters."""
    input:
        #vcf = config["results"] + "{filter_method}/filtered/{dataset}.{mode}.chr{chr}.vcf.gz",
        bcf = config["results"] + "joint_call/biallelic/{dataset}.{mode}.chr{chr}.bcf",
    output:
        bcf = config["results"] + "hard_filtered/pass/{dataset}.{mode}.chr{chr}.bcf",
    params:
        QUAL = CONFIG["QUAL"],
        QD = CONFIG["QD"],
        SOR = CONFIG["SOR"],
        FS = CONFIG["FS"],
        MQ = CONFIG["MQ"],
        MQRankSum = CONFIG["MQRankSum"],
        ReadPosRankSum = CONFIG["ReadPosRankSum"],
    log: config["results"] + "hard_filtered/pass/{dataset}.{mode}.chr{chr}.log",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools view {input.bcf} \
            -e 'QUAL < {params.QUAL} \
                || QD < {params.QD} \
                || SOR < {params.SOR} \
                || FS < {params.FS} \
                || MQ < {params.MQ} \
                || MQRankSum < {params.MQRankSum} \
                || ReadPosRankSum < {params.ReadPosRankSum}' \
            -Ob \
            -o {output.bcf} \
        2> {log}
        """

## --> to phasing.smk

rule merge_filtered_subsets:
    """Merge .vcf file containing SNPs and those containing indels."""
    input:
        bcfs = lambda wildcards: expand(config["results"] + "{filter_method}/pass/{dataset}.{mode}.bcf",
            dataset=wildcards.dataset,
            mode=["SNP", "indel"]),
        tbis = lambda wildcards: expand(config["results"] + "{filter_method}/pass/{dataset}.{mode}.bcf.csi",
            dataset=wildcards.dataset,
            mode=["SNP", "indel"]),
    output:
        config["results"] + "{filter_method}/pass/{dataset}.bcf",
    threads: 2
    resources: nodes = 2
    conda: "../envs/common.yaml"
    shell: """
        bcftools {input.vcfs} \
            --allow-overlaps \
            -Ob \
            -o {output} \
            --threads {threads}"""

# rule create_subset_list:
#     """Create a list of samples for which have the most data. WGS > WES > GBS."""
#     input: all_samples = config["samples"],
#     output: largest_samples = config["results"] + "joint_vcf/largest_samples.list",
#     # sed separates id and sequence type. E.g. WGS12345 -> WGS    12345
#     # awk flips columns
#     # sort rows
#     # merge rows based on first column
#     # take largest id with largest sequence type
#     # sort
#     shell: "sed 's/./&\t/3' {input.all_samples} \
#             | awk -v OFS='\t' '{{print $2,$1}}' \
#             | sort \
#             | awk '$1!=p{{if(p)print s; p=$1; s=$0; next}}{{sub(p,x); s=s $0}} END{{print s}}' \
#             | awk -v OFS='' '{{print $NF,$1}}' \
#             | sort > {output.largest_samples}"        

# rule subset_joint_vcf_largest:
#     """Keep only the sample from each individual with the most data. WGS > WES > GBS. And rename to just animal id.
    
#     To create sample_map from subset_samples:
#     awk -v FS='S' '{print $1"S"$2,$2}' largest_samples.list > largest_samples_map.list
#     """
#     input: vcf = config["results"] + "{filter_method}/pass/{workspace}.{mode}.pass.vcf.gz",  # Need to set up multiallelic rule
#            #subset_samples = config["results"] + "joint_vcf/largest_samples.list",
#            samples_subset = config["resources"] + "samples/largest_samples.list",
#            sample_map = config["resources"] + "samples/largest_samples_map.list",
#     output: subset = config["results"] + "{filter_method}/pass/{workspace}.{mode}.pass.subset.vcf.gz",
#     threads: 1
#     resources: nodes = 1
#     conda: "../envs/common.yaml"
#     # # Ideally reverse order since `reheader` doesn't ahve a -O option to make into .vcf.gz
#     # Shell command does the following:
#     # 1) Subset samples
#     # 2) Rename samples
#     # 3) Keep only variants with an allele count > 0
#     # 4) Combine any remaining multiallelics into a single line.
#     shell: "bcftools view {input.vcf} \
#                 -S {input.samples_subset} \
#                 -Ou \
#             | bcftools reheader \
#                 -s {input.sample_map} \
#             | bcftools view \
#                 --min-ac=1 \
#                 -Ou \
#             | bcftools norm \
#                 -m+any \
#                 -Oz \
#                 -o {output}"

rule subset_joint_vcf:
    """Keep only the sample from each individual with the most data. WGS > WES > GBS. And rename to just animal id.
    
    To create sample_map from subset_samples:
    awk -v FS='S' '{print $1"S"$2,$2}' largest_samples.list > largest_samples_map.list
    """
    input:
        vcf = config["results"] + "{filter_method}/pass/{dataset}.{mode}.pass_only.vcf.gz",  # Need to set up multiallelic rule
        samples_subset = config["resources"] + "samples/largest_samples.list",
        sample_map = config["resources"] + "samples/largest_samples_map.list",
    output:
        subset = config["results"] + "{filter_method}/pass/{dataset}.{mode}.pass_only.subset.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    # # Ideally reverse order since `reheader` doesn't ahve a -O option to make into .vcf.gz
    # Shell command does the following:
    # 1) Subset samples
    # 2) Rename samples
    # 3) Keep only variants with an allele count > 0
    # 4) Combine any remaining multiallelics into a single line.
    shell: """
        bcftools view {input.vcf} \
            -S {input.samples_subset} \
            -Ou \
        | bcftools reheader \
            -s {input.sample_map} \
        | bcftools view \
            --min-ac=1 \
            -Ou \
        | bcftools norm \
            -m+any \
            -Oz \
            -o {output}"""

rule subset_to_GBS_WES_or_WGS:
    """Keep only the sample from each individual with the most data. WGS > WES > GBS. And rename to just animal id.
    
    To create sample_map from subset_samples:
    awk -v FS='S' '{print $1"S"$2,$2}' largest_samples.list > largest_samples_map.list
    """
    #wildcard_constraints: subset = "GBS|WES|WGS"
    input:
        vcf = config["results"] + "{filter_method}/pass/GBS_WES_WGS.{mode}.pass_only.subset.vcf.gz",  # Need to set up multiallelic rule
    output:
        subset = config["results"] + "{filter_method}/pass/{subset}.{mode}.stricter_filter.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    # # Ideally reverse order since `reheader` doesn't ahve a -O option to make into .vcf.gz
    # Shell command does the following:
    # 1) Subset samples
    # 4) Combine any remaining multiallelics into a single line.
    shell: """
        bcftools view {input.vcf} \
            -S {input.samples_subset} \
            -Ou \
        | bcftools norm \
            -m+any \
            -Ov \
        | vcftools \
            --vcf - \
            --max-missing 0.7 \
            -c \
            --recode --recode-INFO-all \
        | bgzip > {output}"""
