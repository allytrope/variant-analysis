# I have a lot of trouble with these rules being run when they aren't supposed to
rule bcf_to_vcfgz:
    """A temporary conversion of a BCF to a gzipped VCF as needed by certain tools."""
    input:
        bcf = "{path}/{name}.bcf"
    output:
        vcfgz = temp("{path}/{name}.vcf.gz"),
        csi = temp("{path}/{name}.vcf.gz.csi")
    conda: "../../envs/common.yaml"
    shell: """
        bcftools view {input} \
            -Oz \
            -o {output.vcfgz}; \
        bcftools index {output.vcfgz} \
        """
ruleorder: bcf_to_vcfgz > csi_index

# rule vcfgz_to_bcf:
#     input: vcfgz = "{path}/{name}.vcf.gz"
#     output: bcf = "{path}/{name}.bcf"
#     conda: "../envs/common.yaml"
#     shell: """
#         bcftools view {input} \
#             -Ob \
#             -o {output.bcf} \
#         """

rule vcf_to_bcf:
    input: vcfgz = "{path}/{name}.vcf"
    output: bcf = "{path}/{name}.bcf"
    params:
        O_arg = "b",
    conda: "../../envs/common.yaml"
    threads: 1
    resources:
        nodes = 1,
        mem_mb = 1000,
    shell: """
        bcftools view {input} \
            -O{params.O_arg} \
            -o {output.bcf} \
        """
# use rule vcf_to_bcf as vcf_to_vcfgz with:
#     output: bcf = temp("{path}/{name}.vcf.gz")
#     params:
#         O_arg = "z",

ruleorder: vcf_to_bcf > bcf_to_vcfgz
