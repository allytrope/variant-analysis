"""Tools for merging structural variants from multiple callers."""

env_path = "../../envs"


rule merge_bcf_samples:
    """Merge individual BCF files."""
    # TODO: Could generalize the suffix
    # wildcard_constraints:
    #     ext = "bcf|vcf|vcf.gz"
    input:
        bcf = expand(config["results"] + "structural_variants/{{tool}}/indivs/{sample}.bcf",
            sample=collect_samples(
                fmt="{batch}/{seq}{indiv}_{library}"
            )
        ),
        csi = expand(config["results"] + "structural_variants/{{tool}}/indivs/{sample}.bcf.csi",
            sample=collect_samples(
                fmt="{batch}/{seq}{indiv}_{library}"
            )
        ),
    output:
        bcf = config["results"] + "structural_variants/{tool}/{dataset}.bcf",
    threads: 1
    resources: nodes = 1
    conda: f"{env_path}/common.yaml"
    shell: """
        bcftools merge {input.bcf} \
            -m none \
            -Ob \
            -o {output.bcf} \
        """

        # | bcftools +fill-tags \
        #     -Ob \
        #     -o {output.bcf} \
        #     -- \
        #     -t AF \

rule sniffles_combined_calling:
    """Use this instead of rule merge_bcf_samples to call SVs from multiple samples at once with sniffles."""
    input:
        snf = expand(config["results"] + "structural_variants/sniffles/indivs/{sample}.snf",
            sample=collect_samples(
                fmt="{batch}/{seq}{indiv}_{library}"
            )
        ),
    output:
        vcf = config["results"] + "structural_variants/sniffles/{dataset}.vcf",
    threads: 1
    resources:
        nodes = 1,
        mem_mb = 14_000,  #20_000
    conda: f"{env_path}/smoove.yaml"
    shell: """
        sniffles \
            --input {input.snf} \
            --vcf {output.vcf} \
            --combine-pctseq 0 \
        """
rule sniffles_to_bcf:
    """Convert sniffles VCF to BCF. This removes any characters including '.' and after it.
    This is necessary because sniffles has been using the file name as the sample name,
    which includes the '.chr{chr}' part."""
    input:
        vcf = config["results"] + "structural_variants/sniffles/{dataset}.vcf",
    output:
        bcf = config["results"] + "structural_variants/sniffles/{dataset}.bcf",
    threads: 1
    resources: nodes = 1
    conda: f"{env_path}/common.yaml"
    shell: """
        bcftools reheader {input.vcf} \
            -N <(paste <(bcftools query -l {input.vcf}) <(bcftools query -l {input.vcf} | cut -d '.' -f 1)) \
        | bcftools view \
            -Ob \
            -o {output.bcf} \
        """
# ruleorder: sniffles_to_bcf > vcf_to_bcf

rule survivor_on_cuteSV:
    """Use SURVIVOR to merge cuteSV calls from multiple samples as was done in Ray et al. paper."""
    input:
        # Can't take BCFs, but can take VCF or gzipped VCF
        vcfs = expand(config["results"] + "structural_variants/cuteSV/indivs/{sample}.vcf.gz",
            sample=collect_samples(
                fmt="{batch}/{seq}{indiv}_{library}"
            )
        ),
    output:
        vcf = config["results"] + "structural_variants/SURVIVOR_from_cuteSV/{dataset}.vcf.gz",
    params:
        # The max_distance and min_length were the only parts specified in the Ray et al. paper
        max_distance = 1000,  # in bases
        min_consensus = 0,
        match_type = 0,  # 0=False, 1=True
        match_strand = 0,  # 0=False, 1=True
        unused = 0,  # Current versions of SURVIVOR do not use this
        min_length = 50,  # Minimum length of SVs
    threads: 1
    resources:
        nodes = 1,
        mem_mb = 30_000,
    conda: f"{env_path}/smoove.yaml"
    shell: """
        SURVIVOR merge <(for FILE in {input.vcfs}; do echo $FILE; done) \
            {params.max_distance} \
            {params.min_consensus} \
            {params.match_type} \
            {params.match_strand} \
            {params.unused} \
            {params.min_length} \
            {output.vcf} \
        """

rule truvari_collapse:
    """Merge variants. Substitute for SURVIVOR."""
    input:
        vcf = config["results"] + "structural_variants/{tool}/vcfgz/{dataset}.vcf.gz",
        # truvari needs requires .tbi specifically, not .csi
        tbi = config["results"] + "structural_variants/{tool}/vcfgz/{dataset}.vcf.gz.tbi",
    output:
        vcf = temp(config["results"] + "structural_variants/truvari_on_{tool}/{dataset}.vcf"),
        collapsed = config["results"] + "structural_variants/truvari_on_{tool}_collapsed/{dataset}.vcf",
    threads: 1
    resources:
        nodes = 1,
        mem_mb = 30_000,
    conda: f"{env_path}/delly2.yaml"
    shell: """
        truvari collapse \
            -i {input.vcf} \
            -o {output.vcf} \
            -c {output.collapsed} \
        """

rule truvari:
    """Comparing SVs using Truvari."""
    input:
        base_calls = config["results"] + "structural_variants/delly/merged/LRS.vcf.gz",  # Must end in .gz and be bgzipped
        base_calls_idx = config["results"] + "structural_variants/delly/merged/LRS.vcf.gz.tbi",
        comp_calls = config["results"] + "structural_variants/delly/merged/4WGS.vcf.gz",  # Must end in .gz and be bgzipped
        comp_calls_idx = config["results"] + "structural_variants/delly/merged/4WGS.vcf.gz.tbi",  # Must end in .gz and be bgzipped
    output:
        config["results"] + "structural_variants/test/test.txt",
        # config["results"] + "structural_variants/truvari/tp-base.vcf.gz",
        # config["results"] + "structural_variants/truvari/tp-comp.vcf.gz",
        # config["results"] + "structural_variants/truvari/fp.vcf.gz",
        # config["results"] + "structural_variants/truvari/fn.vcf.gz",
        # config["results"] + "structural_variants/truvari/summary.json",
        #config["results"] + "structural_variants/truvari/params.json",
    params:
        out_dir = config["results"] + "structural_variants/truvari/",
    threads: 1
    resources: nodes = 1
    conda: f"{env_path}/delly2.yaml"
    shell: """
        truvari bench \
            -b {input.base_calls} \
            -c {input.comp_calls} \
            -o {params.out_dir} \
        """


# ----------------------#
#  SURVIVOR (consensus) #
# ----------------------#

rule SV_consensus:
    """Find consensus between SV callers."""
    input:
        list_of_vcfs = lambda wildcards: expand(config["results"] + "structural_variants/{SV_caller}/merged/{dataset}.vcf",
            SV_caller=["pbsv", "sniffles", "cuteSV"],  # For LRS
            # SV_caller=["delly", "smoove"],  # For WGS (non-LRS)
            dataset=wildcards.dataset),
        # list_of_vcfs = expand(config["results"] + "structural_variants/delly/merged/{dataset}.vcf", dataset=["LRS", "4WGS"]),
    output:
        merged = config["results"] + "structural_variants/SURVIVOR/{dataset}.vcf",
    params:
        max_distance = 1000,  # in bases
        min_consensus = 2,
        match_type = 0,  # 0=False, 1=True
        match_strand = 1,  # 0=False, 1=True
        unused = 0,  # Current versions of SURVIVOR do not use this parameter anymore
        min_length = 0,  # Minimum length of SVs
    threads: 24
    resources: nodes = 24
    conda: f"{env_path}/delly2.yaml"
    shell: """
        SURVIVOR merge <(for FILE in {input.list_of_vcfs}; do echo $FILE; done) \
            {params.max_distance} \
            {params.min_consensus} \
            {params.match_type} \
            {params.match_strand} \
            {params.unused} \
            {params.min_length} \
            {output.merged} \
        """

# rule all_SURVIVOR_combinations:
#     """Summarize overlap between all callers. Requires running every combination of 2 or more callers with SURVIVOR."""
#     input: "",
#     output: "",
#     threads: 24
#     resources: nodes = 24
#     conda: "../envs/delly2.yaml"
#     shell: """

#         """