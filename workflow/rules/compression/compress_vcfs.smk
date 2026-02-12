# I have a lot of trouble with these rules being run when they aren't supposed to
rule bcf_to_vcfgz:
    input: bcf = "{path}/{name}.bcf"
    output: vcfgz = temp("{path}/{name}.vcf.gz")
    priority: -10
    conda: "../../envs/common.yaml"
    shell: """
        bcftools view {input} \
            -Oz \
            -o {output.vcfgz} \
        """

# rule vcfgz_to_bcf:
#     input: vcfgz = "{path}/{name}.vcf.gz"
#     output: bcf = "{path}/{name}.bcf"
#     conda: "../envs/common.yaml"
#     shell: """
#         bcftools view {input} \
#             -Ob \
#             -o {output.bcf} \
#         """