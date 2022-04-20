"""Contain rules for variant recalibration, also called VQSR. Runs after variant_calling.smk as an alternative to hard_filter.smk."""## Recalibration (Alternative to Hard Filtering)---------------
# If variant_calling.smk created a VCF containing WGS data mixed with WES or GBS,
# they should be split for the next step and given different `-an` metrics.

CONFIG = config["variant_calling"]


def modes(mode):
    """Use in rule subset_mode."""
    if mode == "SNP":
        return "-select-type SNP"
    elif mode == "indel":
        return "-select-type INDEL \
                -select-type MIXED"
rule split_by_chromosome_and_mode:
    """Split into SNP- or indel-only .vcf as well as by chromosome (or contigs). The wildcard `mode` can be "SNP" or "indel"."""
    #wildcard_constraints: mode = "SNP|indel",
    wildcard_constraints: mode = "SNP|indel",
    input: split = "{path}/{name}.vcf.gz",
           split_index = "{path}/{name}.vcf.gz.tbi",
    #output: subset = config["results"] + "joint_call/{workspace}.{chrom}.{mode}.vcf.gz",
    output: subset = "{path}/{name}.chr{chr}.{mode}.vcf.gz",
    params: modes = lambda wildcards: modes(wildcards.mode),
    threads: 1
    conda: "../envs/gatk.yaml"
    shell: "gatk SelectVariants \
                -V {input.split} \
                {params.modes} \
                -L {wildcards.chr} \
                -O {output.subset}"

def truth_or_training(string):
    if string == "training":
        return ">="
    elif string == "truth":
        return "<"
rule split_into_truth_and_training:
    """If only one reference .vcf file is available, split into separate truth and training files."""
    wildcard_constraints: resource = "truth|training",
    input: vcf = "{path}/{name}.chr{chr}.{mode}.vcf.gz",
    output: vcf = "{path}/truth_and_training/{name}.chr{chr}.{mode}.{resource}.vcf.gz",
    params: inequality = lambda wildcards: truth_or_training(wildcards.resource),
            threshold = 60,
    conda: "../envs/gatk.yaml"
    shell: """
        gatk VariantFiltration \
            -V {input.vcf} \
            --filter-name "TruAndTrn" \
            --filter-expression "QUAL {params.inequality} {params.threshold}" \
            -O {output.vcf}
        """

rule remove_filtered:
    """Remove variants that have been filtered."""
    input: vcf = "{path}/truth_and_training/{name}.chr{chr}.{mode}.{resource}.vcf.gz",
    output: vcf = "{path}/truth_and_training_removed/{name}.chr{chr}.{mode}.{resource}.vcf.gz",
    conda: "../envs/bio.yaml"
    shell: """
        bcftools view -f PASS {input.vcf} > {output.vcf}
        """

def an_flags(mode):
    """Return `-an` flags and arguments for variant recalibration."""
    if mode == "SNP":
        return "-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP"
    elif mode == "indel":
        return "-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP"
# def resources(mode):
#     """Find corresponding resoures."""
#     if mode == "SNP":
#         pass
#     elif mode == "indel":
#         pass
rule variant_recalibration:
    """Build a recalibration model. Can be run in "SNP" or "indel" mode depending on which string is used for the wildcard `mode`."""
    wildcard_constraints: mode = "SNP|indel",
    input: ref = config["results"] + "ref/ref_genome.fna.gz",
           vcf = config["results"] + "joint_call/renamed/{workspace}.vcf.gz",  # Location to be changed
           truth = lambda wildcards: CONFIG[wildcards.mode]["VQSR_truth_vcf"],
           truth_idx = lambda wildcards: CONFIG[wildcards.mode]["VQSR_truth_vcf"] + ".tbi",
           training = lambda wildcards: CONFIG[wildcards.mode]["VQSR_training_vcf"],
           training_idx = lambda wildcards: CONFIG[wildcards.mode]["VQSR_training_vcf"] + ".tbi",
           #known = lambda wildcards: CONFIG[wildcards.mode]["VQSR_known_vcf"],
           #known_idx = lambda wildcards: CONFIG[wildcards.mode]["VQSR_known_vcf"] + ".tbi",
    output: recal = config["results"] + "VQSR/{workspace}.{mode}.recal",
            recal_idx = config["results"] + "VQSR/{workspace}.{mode}.recal.idx",
            tranches = config["results"] + "VQSR/{workspace}.{mode}.tranches",
            tranches_pdf = config["results"] + "VQSR/{workspace}.{mode}.tranches.pdf",
            fig = config["results"] + "VQSR/{workspace}.{mode}.fig",
            fig_pdf = config["results"] + "VQSR/{workspace}.{mode}.fig.pdf",
    threads: 1
    conda: "../envs/gatk.yaml"
    params: mode = lambda wildcards: wildcards.mode.upper(),
            an_flags = lambda wildcards: an_flags(wildcards.mode),
            truth_prior = lambda wildcards: CONFIG[wildcards.mode]["VQSR_truth_prior"],
            training_prior = lambda wildcards: CONFIG[wildcards.mode]["VQSR_training_prior"],
            #known_prior = lambda wildcards: CONFIG[wildcards.mode]["VQSR_known_prior"],
    # --resource:known,known=true,truth=false,training=false,prior={params.known_prior} {input.known} \
    shell: "gatk --java-options '-Xmx8g' VariantRecalibrator \
                -R {input.ref} \
                -V {input.vcf} \
                -O {output.recal} \
                --resource:truth,known=false,truth=true,training=true,prior={params.truth_prior} {input.truth} \
                --resource:training,known=false,truth=false,training=true,prior={params.training_prior} {input.training} \
                {parameters.an_flags} \
                -mode {params.mode} \
                --tranches-file {output.tranches} \
                --rscript-file {output.fig}"

rule apply_variant_recalibration:
    """Apply recalibration model. Can be run in "SNP" or "indel" mode depending on which string is used for the wildcard `mode`."""
    input: ref = config["results"] + "ref/ref_genome.fna.gz",
           vcf = config["results"] + "joint_call/{workspace}.vcf.gz",
           recal = config["results"] + "VQSR/{workspace}.{mode}.recal",
           tranches = config["results"] + "VQSR/{workspace}.{mode}.tranches",
    output: config["results"] + "joint_call/{workspace}_recalibrated.{mode}.vcf.gz",
    threads: 1
    conda: "../envs/gatk.yaml"
    params: mode = lambda wildcards: wildcards.mode.upper(),
    shell: "gatk --java-options '-Xmx8g' ApplyVQSR \
                -R {input.ref} \
                -V {input.vcf} \
                -O {output} \
                -mode {params.mode} \
                --recal-file {input.recal} \
                --tranches-file {input.tranches}"

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
rule snp_summary:
    """Gives a summary of SNP data from .vcf file."""
    input: vcf = config["results"] + "joint_call/{workspace}_recalibrated.vcf.gz",
           ref = config["results"] + "ref/ref_genome.fna.gz",
    output: txt = config["results"] + "joint_call/{workspace}_snp_summary.txt",
            html = config["results"] + "joint_call/{workspace}_snp_summary.html",
            warnings = config["results"] + "joint_call/{workspace}_snp_summary_warnings.txt",
    params: species = config["ref"]["species"],
    threads: 1
    conda: "../envs/bio.yaml"
    shell: "vep \
    --input_file {input.vcf} \
    --fasta {input.ref} \
    --output_file {output.txt} \
    --stats_text {output.html} \
    --warning_file {output.warnings} \
    --species {params.species} \
    --cache"

