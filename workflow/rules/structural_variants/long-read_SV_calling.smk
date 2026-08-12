### LRS

env_path = "../../envs"

## pbsv

rule pbsv_discover:
    """Reduce PacBio alignments to what is relevant for structural variants using pbsv."""
    input:
        alignment = config["results"] + "alignments/pbmm2/{batch}/{seq}{indiv}_{library}.bam",
        trf = config["resources"] + "trf/rheMac10.trf.bed",
    output:
        svsig = config["results"] + "svsigs/pbsv/indivs/{batch}/{seq}{indiv}_{library}.svsig.gz",
    threads: 16
    resources:
        nodes = 16,
        mem_mb = 8_000,
    conda: f"{env_path}/smoove.yaml"
    #            --region {wildcards.chr} \
    shell: """
        pbsv discover {input.alignment} {output.svsig} \
            --tandem-repeats {input.trf} \
        """

rule pbsv_call:
    """Call structural variants from PacBio using pbsv."""
    # input:
    #     ref = config["ref_fasta"],
    #     svsig = expand(config["results"] + "svsigs/pbsv/{sample}.chr{{chr}}.svsig.gz",
    #         sample=collect_samples(fmt="{batch}/{seq}{indiv}_{library}")),
    #     tbi = expand(config["results"] + "svsigs/pbsv/{sample}.chr{{chr}}.svsig.gz.tbi",
    #         sample=collect_samples(fmt="{batch}/{seq}{indiv}_{library}")),
    # output:
    #     vcf = config["results"] + "structural_variants/pbsv/{dataset}.chr{chr}.vcf.gz",
    input:
        #ref = config["ref_fasta"],
        # TODO: Generalize. This tool won't accept gzipped fasta
        ref = "/master/abagwell/variant-analysis/resources/rhesus/ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa",
        svsig = config["results"] + "svsigs/pbsv/indivs/{batch}/{seq}{indiv}_{library}.svsig.gz",
        tbi = config["results"] + "svsigs/pbsv/indivs/{batch}/{seq}{indiv}_{library}.svsig.gz.tbi",
    output:
        # Marked temporary because this output should be made into a BCF
        vcf = temp(config["results"] + "structural_variants/pbsv/indivs/{batch}/{seq}{indiv}_{library}.vcf"),
    threads: 8
    resources:
        nodes = 8,
        mem_mb = 24_000,  # Tried at 8Gb, but failed
    conda: f"{env_path}/smoove.yaml"
    shell: """
        pbsv call {input.ref} {input.svsig} {output.vcf} \
            --ccs \
            -j {threads}; \
        """

## Sniffles

## Not split by contig
# rule sniffles:
#     """Call structural variants from PacBio using sniffles."""
#     input:
#         ref_fasta = "/master/abagwell/variant-analysis/resources/rhesus/ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa",
#         alignment = config["results"] + "alignments/pbmm2/{batch}/{seq}{indiv}_{library}.bam",
#     output:
#         snf = config["results"] + "structural_variants/sniffles/{batch}/{seq}{indiv}_{library}.snf",
#     threads: 4  # Sniffles uses 4 by default
#     resources:
#         nodes = 4,
#         mem_mb = 20_000,
#     conda: f"{env_path}/smoove.yaml"
# 	shell: """
#         sniffles \
#             --input {input.alignment} \
#             --reference {input.ref_fasta} \
#             --snf {output.snf} \
# 		"""

rule sniffles:
    """Call structural variants from PacBio using sniffles."""
    input:
        ref_fasta = "/master/abagwell/variant-analysis/resources/rhesus/ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa",
        alignment = config["results"] + "alignments/pbmm2/{batch}/{seq}{indiv}_{library}.bam",
    output:
        snf = config["results"] + "structural_variants/sniffles/indivs/{batch}/{seq}{indiv}_{library}.snf",
    threads: 4  # Sniffles uses 4 by default
    resources:
        nodes = 4,
        mem_mb = 14_000,  # This works when I split by chromosome
    conda: f"{env_path}/smoove.yaml"
    #            --contig {wildcards.chr} \
	shell: """
        sniffles \
            --input {input.alignment} \
            --reference {input.ref_fasta} \
            --snf {output.snf} \
		"""

# cuteSV

rule cuteSV:
    input:
        ref_fasta = "/master/abagwell/variant-analysis/resources/rhesus/ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa",
        alignment = config["results"] + "alignments/pbmm2/{batch}/{seq}{indiv}_{library}.bam",
    output:
        # The extra "{seq}{indiv}_{library}" is needed because otherwise, running multiple of this rule
        # will cause errors from each trying to write to the same file
        vcf = temp(config["results"] + "structural_variants/cuteSV/indivs/{batch}/{seq}{indiv}_{library}.vcf"),
    params:
        work_path = config["results"] + "structural_variants/cuteSV/indivs/{batch}/workspace/{seq}{indiv}_{library}",
    threads: 16  # The default
    resources:
        nodes = 16,
        mem_mb = 18_000,
    conda: f"{env_path}/smoove.yaml"
	shell: """
        mkdir {params.work_path} -p; \
        cuteSV {input.alignment} {input.ref_fasta} {output.vcf} {params.work_path} \
            --max_cluster_bias_INS 1000 \
            --diff_ratio_merging_INS 0.9 \
            --max_cluster_bias_DEL 1000 \
            --diff_ratio_merging_DEL 0.5 \
            --genotype \
            -S {wildcards.seq}{wildcards.indiv}_{wildcards.library} \
        """
# SURVIVOR used in between these rules
rule cuteFC:
    """Force-call BAMs using cuteSV VCF that was passed through SURVIVOR."""
    input:
        ref_fasta = "/master/abagwell/variant-analysis/resources/rhesus/ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa",
        alignment = config["results"] + "alignments/pbmm2/{batch}/{seq}{indiv}_{library}.bam",
        #vcf = config["results"] + "structural_variants/SURVIVOR_from_cuteSV/" + config["dataset"] + ".vcf.gz",
        vcf = config["results"] + "structural_variants/truvari_on_cuteSV/" + config["dataset"] + ".vcf.gz",
    output:
        # The extra "{seq}{indiv}_{library}" is needed because otherwise, running multiple of this rule
        # will cause errors from each trying to write to the same file
        vcf = temp(config["results"] + "structural_variants/cuteFC/indivs/{batch}/{seq}{indiv}_{library}.vcf"),    
    params:
        work_path = config["results"] + "structural_variants/cuteFC/indivs/{batch}/workspace/{seq}{indiv}_{library}",
    threads: 16  # The default
    resources:
        nodes = 16,
        mem_mb = 18_000,
    conda: f"{env_path}/smoove.yaml"
    shell: """
        mkdir {params.work_path} -p; \
        cuteFC {input.alignment} {input.ref_fasta} {output.vcf} {params.work_path} \
            -Ivcf {input.vcf} \
            --max_cluster_bias_INS 1000 \
            --diff_ratio_merging_INS 0.9 \d
            --max_cluster_bias_DEL 1000 \
            --diff_ratio_merging_DEL 0.5 \
        """

rule sort_cuteFC:
    """CuteFC doesn't appear to sort by default, so this is needed."""
    input:
        vcf = config["results"] + "structural_variants/cuteFC/indivs/{batch}/{seq}{indiv}_{library}.vcf"
    output:
        bcf = config["results"] + "structural_variants/cuteFC/indivs/{batch}/{seq}{indiv}_{library}.bcf",
    shell: """
        bcftools sort {input} \
            -Ob \
            -o {output.bcf}; \
        """
ruleorder: sort_cuteFC > vcfgz_to_bcf