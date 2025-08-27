"""Rules trimming, alignment, and post-alignment processing of BAMs. Subsequent rules of main workflow are found in variant_calling.smk."""

## Trimming
rule trim_fastp:
    """Trim using fastp."""
    # TODO: Generalize input
    input:
        #reads = lambda wildcards: expand(config["resources"] + "reads/gzipped/{batch}/{seq}{indiv}_{library}_{flowcell_lane}.{read}.fastq.gz",
        reads = lambda wildcards: expand(config["resources"] + "reads/{batch}/{seq}{indiv}_{library}_{flowcell_lane}.{read}.fastq.gz",
            batch=wildcards.batch,
            seq=wildcards.seq,
            indiv=wildcards.indiv,
            library=wildcards.library,
            flowcell_lane=wildcards.flowcell_lane,
            read=["R1", "R2"]),
        #fastq = config["reads"] + "{batch}/{seq}{sample_run}.R1+2.fastq.genozip",
        #fastq = config["results"] + "reads/{batch}/{seq}{sample_run}.R1+R2.fastq",
        ref_genozip = config["compression"]["ref_fasta"],
    output:
        html = config["results"] + "trimmed/html/{batch}/{seq}{indiv}_{library}_{flowcell_lane}.html",
        json = config["results"] + "trimmed/json/{batch}/{seq}{indiv}_{library}_{flowcell_lane}.json",
        trimmed = pipe(config["results"] + "trimmed/reads/{batch}/{seq}{indiv}_{library}_{flowcell_lane}.R1+2.fastq"),
    log:
        config["results"] + "trimmed/log/{batch}/{seq}{indiv}_{library}_{flowcell_lane}.log",
    conda: "../envs/fastp.yaml"
    shell: """
        fastp \
            --in1 {input.reads[0]} \
            --in2 {input.reads[1]} \
            --html {output.html} \
            --json {output.json} \
            --report_title {wildcards.batch}/{wildcards.seq}{wildcards.indiv}_{wildcards.library}_{wildcards.flowcell_lane} \
            --stdout \
            --detect_adapter_for_pe \
            --correction \
        > {output.trimmed} \
        2> {log} \
        """

## Alignment
# rule align_LRS:
#     """Align long read sequencing."""
#     input:
#         mmi = config["ref_fasta"] + ".mmi",
#     output:
#     shell: """
#         minimap2 {input.fastq} \
#             -d {input.mmi} \
#         > {output} \
#         """

rule align:
    """Align .fastq sequence to reference genome.

    Extracts read group information. Then aligns to references genome while adding this read group info.
    Then sorts the output by coordinates."""
    input:
        #fastq_first_line = config["results"] + "reads/{batch}/{sample_run}.first_line.R1+2.fastq",
        ref = config["ref_fasta"],
        ref_indices = multiext(config["ref_fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
        trimmed = config["results"] + "trimmed/reads/{batch}/{seq}{indiv}_{library}_{flowcell_lane}.R1+2.fastq",
    output:
        alignment = pipe(config["results"] + "alignments/raw/{batch}/{seq}{indiv}_{library}_{flowcell_lane}.bam"),
    log:
        config["results"] + "alignments/log/{batch}/{seq}{indiv}_{library}_{flowcell_lane}.log",
    threads: 20
    resources: nodes = 20
    conda: "../envs/align.yaml"
    # First of RG's tags must be SM and last must be PU because of how I have to call the sample names.
    shell: """
        bwa mem {input.ref} {input.trimmed} \
            -p \
            -t {threads} \
            -R '@RG\\tSM:{wildcards.indiv}\\tLB:{wildcards.library}\\tID:{wildcards.indiv}\\tPL:ILLUMINA' \
            2> {log} \
        > {output.alignment} \
        """

rule alignment_postprocessing:
    """Fix mate pairs and mark duplicate reads."""
    wildcard_constraints:
        seq = "WES|WGS|lpWGS",
    input:
        alignment = config["results"] + "alignments/raw/{batch}/{seq}{indiv}_{library}_{flowcell_lane}.bam",
    output:
        alignment = config["results"] + "alignments/markdup/{batch}/{seq}{indiv}_{library}_{flowcell_lane}.bam",
    conda: "../envs/common.yaml"
    threads: 1
    resources: nodes = 1
    # This line can go under each samtools if a specific tmp directory is needed: -T ~/tmp/{rule} \
    shell: """
        samtools sort {input.alignment} \
            -n \
            -u \
        | samtools fixmate - - \
            -m \
            -u \
        | samtools sort - \
            -u \
        | samtools markdup - {output.alignment} \
        """

rule alignment_postprocessing_without_markdup:
    """Fix mate pairs and sort reads. This implementation for AMP and GBS data doesn't mark duplicates
    due to the large false positive rate of duplicates in these sequencing methods."""
    wildcard_constraints:
        seq = "AMP|GBS",
    input:
        alignment = config["results"] + "alignments/raw/{batch}/{seq}{indiv}_{library}_{flowcell_lane}.bam",
    output:
        alignment = config["results"] + "alignments/markdup/{batch}/{seq}{indiv}_{library}_{flowcell_lane}.bam",
    conda: "../envs/common.yaml"
    threads: 1
    resources: nodes = 1
    shell: """
        samtools sort {input.alignment} \
            -n \
            -u \
        | samtools fixmate - - \
            -m \
            -u \
        | samtools sort - \
            -o {output.alignment} \
        """

rule merge_bams:
    """Merge BAMs for same sample.
    This rule is avoided where possible with process substitution,
    but is here for tools that require a seekable BAM."""
    input:
        # bams = collect_runs_from_sample,
        # bams_idx = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_sample(wildcards))),
        bams = lambda wildcards: expand(
            config["results"] + "alignments/markdup/{collect}.bam",
            collect = collect_samples(
                fmt="{batch}/{seq}{indiv}_{library}_{flowcell_lane}",
                col="library",
                val=wildcards.library),
        ),
        bais = lambda wildcards: expand(
            config["results"] + "alignments/markdup/{collect}.bam.bai",
            collect = collect_samples(
                fmt="{batch}/{seq}{indiv}_{library}_{flowcell_lane}",
                col="library",
                val=wildcards.library),
        ),
    output:
        bam = config["results"] + "alignments/merged/{batch}/{seq}{indiv}_{library}.bam",
    conda: "../envs/common.yaml"
    threads: 1
    resources: nodes = 1
    shell: """
        samtools merge {input.bams} \
            -o {output.bam} \
        """