# I have a lot of trouble with these rules being run when they aren't supposed to

rule vcf_to_vcfgz:
    """Compress VCF to a gzipped VCF. Serves as a temporary conversion before conversion to BCF.
    This two step process prevents cyclic dependencies between rules."""
    input: vcf = "{path}/{name}.vcf"
    output:
        vcfgz = temp("{path}/{name}.vcf.gz"),
        csi = temp("{path}/{name}.vcf.gz.csi"),
    conda: "../../envs/common.yaml"
    shell: """
        bcftools view {input} \
            -Oz \
            -o {output.vcfgz}; \
        bcftools index {output.vcfgz} \
        """

# To BCF

rule vcfgz_to_bcf:
    input: vcfgz = "{path}/{name}.vcf.gz"
    output: bcf = "{path}/{name}.bcf"
    conda: "../../envs/common.yaml"
    shell: """
        bcftools view {input} \
            -Ob \
            -o {output.bcf} \
        """

# rule vcf_to_bcf:
#     input: vcf = "{path}/{name}.vcf"
#     output: bcf = "{path}/{name}.bcf"
#     conda: "../envs/common.yaml"
#     shell: """
#         bcftools view {input} \
#             -Ob \
#             -o {output.bcf} \
#         """

# use rule vcf_to_bcf as vcf_to_vcfgz with:
#     output: bcf = temp("{path}/{name}.vcf.gz")
#     params:
#         O_arg = "z",

#ruleorder: vcf_to_bcf > bcf_to_vcfgz

# From BCF

rule bcf_to_vcfgz:
    """A temporary conversion of a BCF to a gzipped VCF as needed by certain tools."""
    input:
        bcf = "{path}/{name}.bcf"
    output:
        vcfgz = temp("{path}/vcfgz/{name}.vcf.gz"),
        csi = temp("{path}/vcfgz/{name}.vcf.gz.csi")
    conda: "../../envs/common.yaml"
    shell: """
        bcftools view {input} \
            -Oz \
            -o {output.vcfgz}; \
        bcftools index {output.vcfgz} \
        """
ruleorder: bcf_to_vcfgz > csi_index

# Because these rules interfere a lot with other rules, they require a lot of rule ordering:
ruleorder: biallelics_by_mode > vcfgz_to_bcf
ruleorder: pass_only_hard_filter > bcf_to_vcfgz
ruleorder: genotype_posteriors > bcf_to_vcfgz
ruleorder: genotype_passing > vcfgz_to_bcf #> vcf_to_bcf

## Genozip

# rule vcf_to_genozip:
#     input:
#         vcf = "{path}/{name}.vcf",
#         ref_genozip = config["compression"]["ref_fasta"],
#     output:
#         genozip = "{path}/{name}.vcf.genozip"
#     conda: "../../envs/common.yaml"
#     threads: 1
#     resources:
#         nodes = 1,
#         mem_mb = 1000,
#     shell: """
#         genozip {input.vcf} \
#             --reference {input.ref_genozip} \
#             --output {output.genozip} \
#         """
# use rule vcf_to_genozip as vcfgz_to_genozip with:
#     input:
#         vcf = "{path}/{name}.vcf.gz",
# use rule vcf_to_genozip as bcf_to_genozip with:
#     input:
#         vcf = "{path}/{name}.bcf",


rule vcfgz_to_genozip:
    input:
        vcfgz = "{path}/{name}.vcf.gz",
        ref_genozip = config["compression"]["ref_fasta"],
    output:
        genozip = "{path}/{name}.vcf.genozip"
    conda: "../../envs/common.yaml"
    threads: 1
    resources:
        nodes = 1,
        mem_mb = 1000,
    shell: """
        genozip {input.vcfgz} \
            --reference {input.ref_genozip} \
            --output {output.genozip} \
        """
# rule bcf_to_genozip:
#     input:
#         bcf = "{path}/{name}.bcf",
#         ref_genozip = config["compression"]["ref_fasta"],
#     output:
#         genozip = "{path}/{name}.vcf.genozip"
#     conda: "../../envs/common.yaml"
#     threads: 1
#     resources:
#         nodes = 1,
#         mem_mb = 1000,
#     shell: """
#         genozip {input.bcf} \
#             --reference {input.ref_genozip} \
#             --output {output.genozip} \
#         """