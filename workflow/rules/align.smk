"""Rules trimming, alignment, and post-alignment processing of BAMs. Subsequent rules of main workflow are found in variant_calling.smk."""

from Bio.Seq import Seq
import pandas as pd

def collect_runs_from_sample(wildcards):
    """Find runs from same sample.
    Can take sample name as wildcards.sample or as wildcards.seq + wildcards.iniv.id."""
    try:
        sample = wildcards.sample
    except AttributeError:
        sample = wildcards.seq + wildcards.indiv_id
    sample_runs = []
    for run in SAMPLE_RUNS:
        # if wildcards.sample == run.split("_")[0]:
        if sample in run:
            sample_runs.append(config["results"] + "alignments/recalibrated/" + run + ".bam")
    return sample_runs

## Trimming
rule cut_adapters_with_i5_i7:
    """Cut out 3' end adapters using i7 and i5 adapter information.
    
    This implementation takes i7 and i5 values from the first line of the .R1.fastq file and assumes that all reads in both R1 and R2 files use those same adapters.
    Additionally, this assumes that the Truseq Dual Index Library was used for the adapters surrounding the i7 and i5 ones.
    These are hard coded under "params". The i7 information is used for finding the adapters in R1 and then i5 for R2."""
    wildcard_constraints:
        seq = "WGS|WES|AMP",
    input:
        reads = lambda wildcards: expand(config["reads"] + "{batch}/{seq}{sample_run}.{read}.fastq.gz",
            batch=wildcards.batch,
            seq=wildcards.seq,
            sample_run=wildcards.sample_run,
            read=["R1", "R2"]),
    output: 
        trimmed = temp(expand(config["results"] + "trimmed/{{batch}}/{{seq}}{{sample_run}}.{read}.fastq.gz",
            read=["R1", "R2"])),
    params: 
        pre_i7 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        post_i7 = "ATCTCGTATGCCGTCTTCTGCTTG",
        pre_i5 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        post_i5 = "GTGTAGATCTCGGTGGTCGCCGTATCATT",
    threads: 4
    resources: nodes = 4
    conda: "../envs/bio.yaml"
    shell: """
        # ADAPTERS=$(gunzip -c {input.reads[0]} | head -n 1 | cut -d " " -f 2 | cut -d ":" -f 4);
        # i7=$(echo $ADAPTERS | cut -d "+" -f 1);
        # i5=$(echo $ADAPTERS | cut -d "+" -f 2);
        R1_END_ADAPTER="{params.pre_i7}";  # ${{i7}}{params.post_i7}";
        R2_END_ADAPTER="{params.pre_i5}";  # ${{i5}}{params.post_i5}";
        cutadapt {input.reads} \
            -a $R1_END_ADAPTER \
            -A $R2_END_ADAPTER \
            --cores {threads} \
            -o {output.trimmed[0]} \
            -p {output.trimmed[1]}"""

rule cut_adapters_with_i5_i7_genozipped:
    """Cut out 3' end adapters using i7 and i5 adapter information from genozipped FASTQs.
    
    This implementation takes i7 and i5 values from the first line of the .R1.fastq file and assumes that all reads in both R1 and R2 files use those same adapters.
    Additionally, this assumes that the Truseq Dual Index Library was used for the adapters surrounding the i7 and i5 ones.
    These are hard coded under "params". The i7 information is used for finding the adapters in R1 and then i5 for R2."""
    wildcard_constraints:
        seq = "WGS|WES|AMP",
    input:
        reads = config["reads"] + "{batch}/{seq}{sample_run}.R1+2.fastq.genozip",
        genozip_ref = config["compression"]["ref_fasta"],
    output: 
        trimmed = temp(expand(config["results"] + "trimmed/{{batch}}/{{seq}}{{sample_run}}.{read}.fastq.gz",
            read=["R1", "R2"])),
    params: 
        pre_i7 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        post_i7 = "ATCTCGTATGCCGTCTTCTGCTTG",
        pre_i5 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        post_i5 = "GTGTAGATCTCGGTGGTCGCCGTATCATT",
    threads: 4
    resources: nodes = 4
    conda: "../envs/bio.yaml"
    shell: """
        R1_END_ADAPTER="{params.pre_i7}";  # ${{i7}}{params.post_i7}";
        R2_END_ADAPTER="{params.pre_i5}";  # ${{i5}}{params.post_i5}";
        cutadapt <(genocat {input.reads} --reference {input.genozip_ref} --R1) <(genocat {input.reads} --reference {input.genozip_ref} --R2) \
            -a $R1_END_ADAPTER \
            -A $R2_END_ADAPTER \
            --cores {threads} \
            -o {output.trimmed[0]} \
            -p {output.trimmed[1]}"""

# rule cut_adapters_with_i5_i7_AMP:
#     """Cut out 3' end adapters using i7 and i5 adapter information.
    
#     This implementation takes i7 and i5 values from a barcodes file.
#     Additionally, this assumes that the Truseq Dual Index Library was used for the adapters surrounding the i7 and i5 ones.
#     These are hard coded under "params". The i7 information is used for finding the adapters in R1 and then i5 for R2."""
#     wildcard_constraints:
#         seq = "AMP",
#     input: 
#         reads = lambda wildcards: expand(config["reads"] + "{seq}{organism_id}.{read}.fastq.gz",
#             seq=wildcards.seq,
#             organism_id=wildcards.organism_id,
#             read=["R1", "R2"]),
#         barcodes = config["barcodes"]["AMP_barcodes"],
#     output: 
#         trimmed = expand(config["results"] + "trimmed/{{seq}}{{organism_id}}.{read}.fastq.gz",
#             read=["R1", "R2"]),
#     params: 
#         pre_i7 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
#         post_i7 = "ATCTCGTATGCCGTCTTCTGCTTG",
#         pre_i5 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
#         post_i5 = "GTGTAGATCTCGGTGGTCGCCGTATCATT",
#     threads: 4
#     resources: nodes = 4
#     conda: "../envs/bio.yaml"
#     shell: """
#         i7=$(grep -P "^{wildcards.organism_id}\t" {input.barcodes} | cut -f 2);
#         i5=$(grep -P "^{wildcards.organism_id}\t" {input.barcodes} | cut -f 3);
#         R1_END_ADAPTER="{params.pre_i7}${{i7}}{params.post_i7}";
#         R2_END_ADAPTER="{params.pre_i5}${{i5}}{params.post_i5}";
#         cutadapt {input.reads} \
#             -a $R1_END_ADAPTER \
#             -A $R2_END_ADAPTER \
#             --cores {threads} \
#             -o {output.trimmed[0]} \
#             -p {output.trimmed[1]}"""

## This section started to error. Probably due to lack of barcodes file.
# def find_barcode(wildcards, input):
#     """Pull sample barcode from barcodes file, which contains reverse complements for R2."""
#     df = pd.read_table(input.barcodes, header=None, names=["organism_id", "barcodes"])
#     barcode = df.set_index("organism_id")["barcodes"].to_dict()[wildcards.sample_run]
#     return str(Seq(barcode).reverse_complement())
# rule cut_GBS_adapters:
#     """Cut out 3' end adapters from GBS reads.

#     This method doesn't rely on i7 or i5 values given in the header lines. Instead, R1 and R2 reads have the adapters given under params.
#     After both of these adapters, there is a pseudoconsensus sequence, but left out since Cutadapt will remove anything after the given adapters."""
#     wildcard_constraints:
#         seq = "GBS",
#     input:
#         reads = lambda wildcards: expand(config["reads"] + "{seq}{sample_run}.{read}.fastq.gz",
#             seq=wildcards.seq,
#             sample_run=wildcards.sample_run,
#             read=["R1", "R2"]),
#         barcodes = config["barcodes"]["GBS_barcodes"],
#     output: 
#         trimmed = expand(config["results"] + "trimmed/{{seq}}{{sample_run}}.{read}.fastq.gz",
#             read=["R1", "R2"]),
#     params: 
#         barcode = lambda wildcards, input: find_barcode(wildcards, input),
#         R1_adapter = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",
#         R2_adapter = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
#     threads: 4
#     resources: nodes = 4
#     conda: "../envs/bio.yaml"
#     shell: """
#         ADAPTERS=$(gunzip -c {input.reads[0]} | head -n 1 | cut -d " " -f 2 | cut -d ":" -f 4);
#         R1_END_ADAPTER="{params.R1_adapter}";
#         R2_END_ADAPTER="{params.barcode}{params.R2_adapter}";
#         cutadapt {input.reads} \
#             -a $R1_END_ADAPTER \
#             -A $R2_END_ADAPTER \
#             --cores {threads} \
#             -o {output.trimmed[0]} \
#             -p {output.trimmed[1]}"""

## Alignment

rule align_LRS:
    """Align long read sequencing."""
    input:
        mmi = config["ref_fasta"] + ".mmi",
    output:
    shell: """
        minimap2 {input.fastq} \
            -d {input.mmi} \
        > {output} \
        """

rule align:
    """Align .fastq sequence to reference genome.

    Extracts read group information. Then aligns to references genome while adding this read group info.
    Then sorts the output by coordinates."""
    input:
        ref = config["ref_fasta"],
        idxs = multiext(config["ref_fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
        trimmed = lambda wildcards: expand(config["results"] + "trimmed/{batch}/{sample_run}.{read}.fastq.gz",
            batch=wildcards.batch,
            sample_run=wildcards.sample_run,
            read=["R1", "R2"]),
    output:
        alignment = temp(config["results"] + "alignments/raw/{batch}/{sample_run}.bam"),
    threads: 24
    resources: nodes = 24
    conda: "../envs/bio.yaml"
    # First of RG's tags must be SM and last must be PU because of how I have to call the sample names.
    shell: """
        bash workflow/scripts/align.sh {input.trimmed[0]} {input.trimmed[1]} {input.ref} {output} {threads} {wildcards.sample_run}"""

rule alignment_postprocessing:
    """Fix mate pairs and mark duplicate reads."""
    wildcard_constraints:
        seq = "WES|WGS",
    input:
        alignment = config["results"] + "alignments/raw/{batch}/{seq}{sample_run}.bam",
    output:
        alignment = temp(config["results"] + "alignments/markdup/{batch}/{seq}{sample_run}.bam"),
    conda: "../envs/bio.yaml"
    threads: 1
    resources: nodes = 1
    # This line can go under each samtools            -T ~/tmp/{rule} \
    shell: """
        samtools sort {input.alignment} \
            -n \
        | samtools fixmate - - \
            -m \
        | samtools sort - \
        | samtools markdup - {output.alignment} \
        """

rule alignment_postprocessing_without_markdup:
    """Fix mate pairs and sort reads. This implementation for AMP and GBS data doesn't mark duplicates
    due to the large false positive rate of duplicates in these sequencing methods."""
    wildcard_constraints:
        seq = "AMP|GBS",
    input:
        alignment = config["results"] + "alignments/raw/{batch}/{seq}{sample_run}.bam",
    output:
        alignment = temp(config["results"] + "alignments/markdup/{batch}/{seq}{sample_run}.bam"),
    conda: "../envs/bio.yaml"
    threads: 1
    resources: nodes = 1
    shell: """
        samtools sort {input.alignment} \
            -n \
        | samtools fixmate - - \
            -m \
        | samtools sort - \
            -o {output.alignment} \
        """

## Base recalibration
rule base_recalibration:
    """Create recalibration table for correcting systemic error in base quality scores."""
    input: 
        ref = config["ref_fasta"],
        ref_idx = config["ref_fasta"] + ".fai",
        ref_gzi = config["ref_fasta"] + ".gzi",
        ref_dict = ".".join(config["ref_fasta"].split(".")[:-2]) + ".dict",
        bam = config["results"] + "alignments/markdup/{batch}/{sample_run}.bam",
        known_variants = config["BQSR_known_variants"],
        indexed_known_variants = config["BQSR_known_variants"] + ".tbi",
    output:
        config["results"] + "alignments/recal_tables/{batch}/{sample_run}.BQSR.recal",
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    shell: """
        gatk --java-options '-Xmx8g' BaseRecalibrator \
            -R {input.ref} \
            -I {input.bam} \
            --known-sites {input.known_variants} \
            -O {output} \
        """

rule apply_base_recalibration:
    """Correct systemic error in base quality scores."""
    input: 
        ref = config["ref_fasta"],
        ref_idx = config["ref_fasta"] + ".fai",
        ref_dict = ".".join(config["ref_fasta"].split(".")[:-2]) + ".dict",
        bam = config["results"] + "alignments/markdup/{batch}/{sample_run}.bam",
        recal = config["results"] + "alignments/recal_tables/{batch}/{sample_run}.BQSR.recal",
    output:
        config["results"] + "alignments/recalibrated/{batch}/{sample_run}.bam",
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    shell: """
        gatk --java-options '-Xmx8g' ApplyBQSR \
            -R {input.ref} \
            -I {input.bam} \
            --bqsr-recal-file {input.recal} \
            -O {output} \
            --create-output-bam-index false; \
        """

rule merge_bams:
    """Merge BAMs for same sample.
    This rule is avoided where possible with process substitution,
    but for tools that require a seekable BAM."""
    input:
        bams = collect_runs_from_sample,
        bams_idx = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_sample(wildcards))),
    output:
        #bam = temp(config["results"] + "alignments/merged/{batch}/{seq}{sample_run}.bam"),
        bam = config["results"] + "alignments/merged/{sample}.bam",  # Just for `rhesus` project until I update it
    conda: "../envs/bio.yaml"
    threads: 1
    resources: nodes = 1
    shell: """
        samtools merge {input.bams} \
            -o {output.bam} \
        """
