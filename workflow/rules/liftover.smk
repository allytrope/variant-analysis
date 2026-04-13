"""Rules related to lifting over variants from one reference genome to another."""

# TODO: Finish this
# rule bcftools_liftover:
#     """Liftover using BCFtools/liftover."""
#     input:
#         bcf = '',
#         src_fasta = '',
#         dst_fasta = '/master/abagwell/variant-analysis/resources/rhesus/annotations/hg38.fa.gz',
#         chain = '',
#     output:
#         bcf = '',
#     shell: """
#         bcftools +liftover {input.bcf} \
#             --src-fasta-ref {input.src_fasta} \
#             --fasta-ref {input.dst_fasta} \
#             --chain {inpu.chain} \
#             -Ob \
#             -o {output.bcf} \
#         """

rule liftover_VCF:
    """Liftover a VCF to a new reference genome."""
    #rheMac10ToHg38
    input:
        #bcf = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.bcf",
        #csi = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.bcf.csi",
        vcf = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz.tbi",
        # vcf = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.bcf",
        # tbi = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.bcf.csi",
        #TODO: Generalize
        chain = "/master/abagwell/variant-analysis/resources/rhesus/annotations/{liftover}.over.chain",
        #chain = "/master/abagwell/variant-analysis/resources/rhesus/annotations/{liftover}.over.chain",
        #ref_fasta = config["ref_fasta"],
        #ref_idx = config["ref_fasta"] + ".fai",
        #ref_dict = ".".join(config["ref_fasta"].split(".")[:-2]) + ".dict",
        target_fasta = "/master/abagwell/variant-analysis/resources/rhesus/annotations/hg38.fa.gz",
        target_fai = "/master/abagwell/variant-analysis/resources/rhesus/annotations/hg38.fa.gz.fai",
        target_dict = "/master/abagwell/variant-analysis/resources/rhesus/annotations/hg38.dict",
    output:
        vcf = config["results"] + "liftover/{liftover}/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz",
        reject = config["results"] + "liftover/{liftover}/rejected/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz",
    priority: 1
    threads: 1
    resources:
        nodes = 1,
        mem_mb = 360_000 #64_000,  # 32_000 still caused memory errors
    conda: "../envs/gatk.yaml"
    shell: """
        gatk --java-options "-Xmx128G" LiftoverVcf \
            -I {input.vcf} \
            -O {output.vcf} \
            --REJECT {output.reject} \
            --CHAIN {input.chain} \
            -R {input.target_fasta} \
        """