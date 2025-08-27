"""Contain rules for variant recalibration, also called VQSR. Runs after variant_calling.smk as an alternative to hard_filter.smk."""## Recalibration (Alternative to Hard Filtering)---------------
# If variant_calling.smk created a VCF containing WGS data mixed with WES or GBS,
# they should be split for the next step and given different `-an` metrics.


## Prepare truth and training from single .vcf (for when separate files don't exist)
def truth_or_training(string):
    if string == "training":
        return ">="
    elif string == "truth":
        return "<"
rule split_into_truth_and_training:
    """If only one reference .vcf file is available, split into separate truth and training files."""
    wildcard_constraints: resource = "truth|training",
    input: vcf = "{path}/{name}.vcf.gz",
    output: vcf = "{path}/truth_and_training/{name}.{resource}.vcf.gz",
    params: inequality = lambda wildcards: truth_or_training(wildcards.resource),
            threshold = 60,
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    shell: """
        gatk VariantFiltration \
            -V {input.vcf} \
            --filter-name "TruAndTrn" \
            --filter-expression "QUAL {params.inequality} {params.threshold}" \
            -O {output.vcf}
        """

rule remove_filtered_from_truth_or_training:
    """Remove variants that have been filtered."""
    input: vcf = "{path}/truth_and_training/{name}.chr{chr}.{mode}.{resource}.vcf.gz",
    output: vcf = "{path}/truth_and_training_removed/{name}.chr{chr}.{mode}.{resource}.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools view -f PASS {input.vcf} -Oz > {output.vcf}
        """

## Variant recalibration
def an_flags(mode):
    """Return `-an` flags and arguments for variant recalibration."""
    if mode == "SNP":
        return "-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP"
    elif mode == "indel":
        return "-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP"
rule variant_recalibration:
    """Build a recalibration model. Can be run in "SNP" or "indel" mode depending on which string is used for the wildcard `mode`.
    Note that this rule requires more memory than most."""
    input:
        ref = config["ref_fasta"],
        ref_idx = config["ref_fasta"] + ".fai",
        ref_dict = ".".join(config["ref_fasta"].split(".")[:-2]) + ".dict",  # Replaces the ".fna.gz" ending with ".dict"
        vcf = config["results"] + "joint_call/split/{dataset}.{mode}.biallelic.vcf.gz",  # Location to be changed
        vcf_idx = config["results"] + "joint_call/split/{dataset}.{mode}.biallelic.vcf.gz.tbi",
        truth = lambda wildcards: config["variant_calling"][wildcards.mode]["VQSR_truth_vcf"],
        truth_idx = lambda wildcards: config["variant_calling"][wildcards.mode]["VQSR_truth_vcf"] + ".tbi",
        training = lambda wildcards: config["variant_calling"][wildcards.mode]["VQSR_training_vcf"],
        training_idx = lambda wildcards: config["variant_calling"][wildcards.mode]["VQSR_training_vcf"] + ".tbi",
        #known = lambda wildcards: CONFIG[wildcards.mode]["VQSR_known_vcf"],
        #known_idx = lambda wildcards: CONFIG[wildcards.mode]["VQSR_known_vcf"] + ".tbi",
    output:
        recal = config["results"] + "VQSR/model/{dataset}.{mode}.biallelic.recal",
        recal_idx = config["results"] + "VQSR/model/{dataset}.{mode}.biallelic.recal.idx",
        tranches = config["results"] + "VQSR/model/{dataset}.{mode}.biallelic.tranches",
        tranches_pdf = config["results"] + "VQSR/model/{dataset}.{mode}.biallelic.tranches.pdf",
        fig = config["results"] + "VQSR/model/{dataset}.{mode}.biallelic.fig",
        fig_pdf = config["results"] + "VQSR/model/{dataset}.{mode}.biallelic.fig.pdf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    params:
        mode = lambda wildcards: wildcards.mode.upper(),
        an_flags = lambda wildcards: an_flags(wildcards.mode),
        truth_prior = lambda wildcards: config["variant_calling"][wildcards.mode]["VQSR_truth_prior"],
        training_prior = lambda wildcards: config["variant_calling"][wildcards.mode]["VQSR_training_prior"],
        #known_prior = lambda wildcards: CONFIG[wildcards.mode]["VQSR_known_prior"],
        tmp_dir = config["results"] + "VQSR/model/tmp/",
    # --resource:known,known=true,truth=false,training=false,prior={params.known_prior} {input.known} \
    shell: """
        gatk --java-options '-Xmx24g' VariantRecalibrator \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.recal} \
            --resource:truth,known=false,truth=true,training=true,prior={params.truth_prior} {input.truth} \
            --resource:training,known=false,truth=false,training=true,prior={params.training_prior} {input.training} \
            {params.an_flags} \
            -mode {params.mode} \
            --tranches-file {output.tranches} \
            --rscript-file {output.fig} \
            --tmp-dir ~/tmp/{rule}"""

rule apply_variant_recalibration:
    """Apply recalibration model. Can be run in "SNP" or "indel" mode depending on which string is used for the wildcard `mode`."""
    input:
        ref = config["ref_fasta"],
        ref_idx = config["ref_fasta"] + ".fai",
        vcf = config["results"] + "joint_call/split/{dataset}.{mode}.biallelic.vcf.gz",
        vcf_idx = config["results"] + "joint_call/split/{dataset}.{mode}.biallelic.vcf.gz.tbi",
        recal = config["results"] + "VQSR/model/{dataset}.{mode}.biallelic.recal",
        tranches = config["results"] + "VQSR/model/{dataset}.{mode}.biallelic.tranches",
    output:
        config["results"] + "VQSR/filter_applied/{dataset}.{mode}.biallelic.recalibrated.vcf.gz",
    params:
        mode = lambda wildcards: wildcards.mode.upper(),
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    shell: """
        gatk --java-options '-Xmx24g' ApplyVQSR \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output} \
            -mode {params.mode} \
            --recal-file {input.recal} \
            --tranches-file {input.tranches}"""

rule pass_only:
    """Remove variants that have been filtered."""
    wildcard_constraints: mode = "SNP|indel",
    input:
        vcf = config["results"] + "{filter_method}/filter_applied/{dataset}.{mode}.biallelic.recalibrated.vcf.gz",
    output:
        vcf = config["results"] + "{filter_method}/pass_only/{dataset}.{mode}.biallelic.recalibrated.pass_only.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools view {input.vcf} \
            -f PASS \
            -Oz \
            -o {output.vcf}"""


# Under construction
# rule merge_SNPs_and_indels:
#     input: SNPs_vcf = config["results"] + "joint_call/{workspace}_recalibrated.SNP.vcf.gz",
#            indels_vcf = config["results"] + "joint_call/{workspace}_recalibrated.indel.vcf.gz",
#     output: "",
#     shell: ""

# CHANGE TO WORK FOR MMUL_10 DATA
'''
rule snp_summary:
    """Gives a summary of SNP data from .vcf file."""
    input: vcf=config["results"] + "joint_call/recalibrated_joint_call.vcf.gz"
    output: config["results"] + "joint_call/snp_summary.txt"
    shell: "vep \
    --custom /home/flow/ensembl-vep/Mmul_8.0.1.92.chr.gff3.gz,,gff \
    --fasta Mmul_8.0.1.chromosome.fa \
    --gff /home/flow/ensembl-vep/Mmul_8.0.1.92.chr.gff3.gz \
    --input_file {input.vcf} \
    --output_file {output} \
    --stats_text "
'''

# Download cache for species of interest first
# rule snp_summary:
#     """Gives a summary of SNP data from .vcf file."""
#     input: vcf = config["results"] + "joint_call/{workspace}_recalibrated.vcf.gz",
#            ref = config["results"] + "ref/ref_genome.fna.gz",
#     output: txt = config["results"] + "joint_call/{workspace}_snp_summary.txt",
#             html = config["results"] + "joint_call/{workspace}_snp_summary.html",
#             warnings = config["results"] + "joint_call/{workspace}_snp_summary_warnings.txt",
#     params: species = config["ref"]["species"],
#     threads: 1
#     conda: "../envs/common.yaml"
#     shell: "vep \
#     --input_file {input.vcf} \
#     --fasta {input.ref} \
#     --output_file {output.txt} \
#     --stats_text {output.html} \
#     --warning_file {output.warnings} \
#     --species {params.species} \
#     --cache"

