"""Contain rules for hard filtering. Runs after variant_calling.smk and as an alternative to variant_recalibration.smk."""

CONFIG = config["hard_filter"]

def filters(mode):
    """Use in rule hard_filter."""
    if mode == "SNP":
        return f"""-filter "QUAL < {CONFIG["QUAL"]}" --filter-name "QUAL{CONFIG["QUAL"]}" \
                  -filter "QD < {CONFIG["QD"]}" --filter-name "QD{CONFIG["QD"]}" \
                  -filter "SOR > {CONFIG["SOR"]}" --filter-name "SOR{CONFIG["SOR"]}" \
                  -filter "FS > {CONFIG["FS"]}" --filter-name "FS{CONFIG["FS"]}" \
                  -filter "MQ < {CONFIG["MQ"]}" --filter-name "MQ{CONFIG["MQ"]}" \
                  -filter "MQRankSum < {CONFIG["MQRankSum"]}" --filter-name "MQRankSum{CONFIG["MQRankSum"]}" \
                  -filter "ReadPosRankSum < {CONFIG["ReadPosRankSum"]}" --filter-name "ReadPosRankSum{CONFIG["ReadPosRankSum"]}" """
    elif mode == "indel":
        return """-filter "QUAL < 30.0" --filter-name "QUAL30" \
                  -filter "QD < 2.0" --filter-name "QD2" \
                  -filter "FS > 200.0" --filter-name "FS200" \
                  -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" """
rule hard_filter:
    """Hard filter .vcf file."""
    input:
        vcf = config["results"] + "joint_call/biallelic/{dataset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "joint_call/biallelic/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
    output:
        filtered = config["results"] + "hard_filtered/filtered/{dataset}.{mode}.chr{chr}.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    params:
        filters = lambda wildcards: filters(wildcards.mode),
    shell: """
        gatk --java-options '-Xmx8g' VariantFiltration \
            -V {input.vcf} \
            -O {output.filtered} \
            {params.filters} """

rule pass_only_hard_filter:
    """Remove variants that have been filtered."""
    input:
        vcf = config["results"] + "{filter_method}/filtered/{dataset}.{mode}.chr{chr}.vcf.gz",
    output:
        vcf = config["results"] + "{filter_method}/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools view {input.vcf} \
            -f PASS \
            -Oz \
            -o {output.vcf}"""

## --> to phasing.smk

rule merge_filtered_subsets:
    """Merge .vcf file containing SNPs and those containing indels."""
    input:
        vcfs = lambda wildcards: expand(config["results"] + "{filter_method}/pass/{dataset}.{mode}.vcf.gz",
            dataset=wildcards.dataset,
            mode=["SNP", "indel"]),
        tbis = lambda wildcards: expand(config["results"] + "{filter_method}/pass/{dataset}.{mode}.vcf.gz.tbi",
            dataset=wildcards.dataset,
            mode=["SNP", "indel"]),
    output:
        config["results"] + "{filter_method}/pass/{dataset}.vcf.gz",
    threads: 2
    resources: nodes = 2
    conda: "../envs/bio.yaml"
    shell: """
        bcftools  {input.vcfs} \
            --allow-overlaps \
            -Oz \
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
#     conda: "../envs/bio.yaml"
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
    conda: "../envs/bio.yaml"
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
    conda: "../envs/bio.yaml"
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
