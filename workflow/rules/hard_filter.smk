"""Contain rules for hard filtering. Runs after variant_calling.smk and as an alternative to variant_recalibration.smk."""


def filters(mode):
    """Use in rule hard_filter."""
    if mode == "SNP":
        return """-filter "QUAL < 30.0" --filter-name "QUAL30" \
                  -filter "QD < 2.0" --filter-name "QD2" \
                  -filter "SOR > 3.0" --filter-name "SOR3" \
                  -filter "FS > 60.0" --filter-name "FS60" \
                  -filter "MQ < 40.0" --filter-name "MQ40" \
                  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
                  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" """
    elif mode == "indel":
        return """-filter "QUAL < 30.0" --filter-name "QUAL30" \
                  -filter "QD < 2.0" --filter-name "QD2" \
                  -filter "FS > 200.0" --filter-name "FS200" \
                  -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" """
rule hard_filter:
    """Hard filter .vcf file."""
    wildcard_constraints: mode = "SNP|indel",
    input: vcf = config["results"] + "{path}/{workspace}.{mode}.vcf.gz",
           vcf_idx = config["results"] + "{path}/{workspace}.{mode}.vcf.gz.tbi",
    output: filtered = config["results"] + "{path}/{workspace}.{mode}.filtered.vcf.gz",
    threads: 1
    conda: "../envs/gatk.yaml"
    params: filters = lambda wildcards: filters(wildcards.mode),
    shell: """gatk --java-options '-Xmx8g' VariantFiltration \
                -V {input.vcf} \
                -O {output.filtered} \
                {params.filters} """

### These work for both hard filtering and VQSR

rule pass_only_hard_filter:
    """Remove variants that have been filtered."""
    wildcard_constraints: mode = "SNP|indel",
    input: vcf = config["results"] + "{filter_method}/filter_applied/{workspace}.{mode}.filtered.vcf.gz",
    output: vcf = config["results"] + "{filter_method}/pass_only/{workspace}.{mode}.pass_only.vcf.gz",
    conda: "../envs/bio.yaml"
    shell: """
        bcftools view {input.vcf} \
            -f PASS \
            -Oz \
            -o {output.vcf}"""

rule merge_filtered_subsets:
    """Merge .vcf file containing SNPs and those containing indels."""
    input: vcfs = lambda wildcards: expand(config["results"] + "{filter_method}/pass_only/{workspace}.{mode}.pass_only.vcf.gz", workspace=wildcards.workspace, mode=["SNP", "indel"]),
           tbi_vcfs = lambda wildcards: expand(config["results"] + "{filter_method}/pass_only/{workspace}.{mode}.pass_only.vcf.gz.tbi", workspace=wildcards.workspace, mode=["SNP", "indel"]),
    output: config["results"] + "{filter_method}/pass_only/{workspace}.pass_only.vcf.gz",
    threads: 2
    conda: "../envs/bio.yaml"
    shell: "bcftools concat {input.vcfs} \
                --allow-overlaps \
                -Oz \
                -o {output} \
                --threads {threads}"

# ------

# rule chromosome_remapping:
#     """Change names of contigs in .vcf according to given mapping.
#     For example, can be used to change project contig ids to those of RefSeq.
#     Or to change be between form {1, 2, 3,...} and {chr1, chr2, chr3,...}."""
#     input: vcf = config["results"] + "joint_call/{workspace}.filtered.vcf.gz",
#            map = CONFIG["chromosome_remapping"],
#     output: config["results"] + "joint_call/{workspace}.filtered.remap.vcf.gz",
#     threads: 1
#     #conda: "../envs/bio.yaml"
#     shell: "bcftools annotate {input.vcf} \
#                 --rename-chrs {input.map} \
#                 -Ou \
#                 | bcftools sort \
#                 -Oz \
#                 -o {output}"

# rule create_subset_list:
#     """Create a list of samples for which have the most data. WGS > WES > GBS."""
#     input: all_samples = config["resources"] + "samples/samples_with_type.list", #  A list of sample names such as "WGS12235" or "GBS14565". Not including those with underscores.
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
#             | awk -v OFS='' '{print $NF,$1}' \
#             | sort > {output.largest_samples}"        

rule subset_joint_vcf:
    """Keep only the sample from each individual with the most data. WGS > WES > GBS. And rename to just animal id.
    
    To create sample_map from subset_samples:
    awk -v FS='S' '{print $1"S"$2,$2}' largest_samples.list > largest_samples_map.list
    """
    input: vcf = config["results"] + "{filter_method}/pass_only/{workspace}.{mode}.pass_only.vcf.gz",  # Need to set up multiallelic rule
           #subset_samples = config["results"] + "joint_vcf/largest_samples.list",
           samples_subset = config["resources"] + "samples/largest_samples.list",
           sample_map = config["resources"] + "samples/largest_samples_map.list",
    output: subset = config["results"] + "{filter_method}/pass_only/{workspace}.{mode}.pass_only.subset.vcf.gz",
    # # Ideally reverse order since `reheader` doesn't ahve a -O option to make into .vcf.gz
    # Shell command does the following:
    # 1) Subset samples
    # 2) Rename samples
    # 3) Keep only variants with an allele count > 0
    # 4) Combine any remaining multiallelics into a single line.
    shell: "bcftools view {input.vcf} \
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
                -o {output}"