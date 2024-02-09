"""Contain workflow for performing variant analysis from paired-end FASTQ files."""

from Bio.Seq import Seq
#import gzip
import os
import pandas as pd


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
    shell: """
        samtools sort {input.alignment} \
            -n \
            -T ~/tmp/{rule} \
        | samtools fixmate - - \
            -m \
        | samtools sort - \
            -T ~/tmp/{rule} \
        | samtools markdup - {output.alignment} \
            -T ~/tmp/{rule} \
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
            -T ~/tmp/{rule} \
        | samtools fixmate - - \
            -m \
        | samtools sort - \
            -T ~/tmp/{rule} \
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
    # Note tmp directory must already exist
    shell: """
        gatk --java-options '-Xmx8g' BaseRecalibrator \
            -R {input.ref} \
            -I {input.bam} \
            --known-sites {input.known_variants} \
            -O {output} \
            --tmp-dir ~/tmp/{rule}"""

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
            --tmp-dir ~/tmp/{rule} \
            --create-output-bam-index false; \
        """

# def collect_runs_from_sample(wildcards):
#     """Find runs from same sample."""
#     sample_runs = []
#     for run in SAMPLE_RUNS:
#         if wildcards.sample in run:
#             sample_runs.append(config["results"] + "alignments/recalibrated/" + run + ".bam")
#     return sample_runs
# rule merge_sample_runs:
#     """Merge runs from same sample."""
#     input:
#         bams = collect_runs_from_sample
#     output:
#         merged = config["results"] + "alignments/merged/{sample}.bam",
#     threads: 1
#     resources: nodes = 1
#     conda: "../envs/bio.yaml"
#     shell: """
#         samtools merge {input.bams} \
#             -o {output.merged} \
#         """

## Variant calling
def collect_runs_from_sample(wildcards):
    """Find runs from same sample."""
    sample_runs = []
    for run in SAMPLE_RUNS:
        # if wildcards.sample == run.split("_")[0]:
        if wildcards.sample_run in run:
            sample_runs.append(config["results"] + "alignments/recalibrated/" + run + ".bam")
    return sample_runs
rule call_variants:
    """Call variants to make VCF file."""
    input:
        ref = config["ref_fasta"],
        fai = config["ref_fasta"] + ".fai",
        dict = ".".join(config["ref_fasta"].split(".")[:-2]) + ".dict",
        #bam = config["results"] + "alignments/merged/{sample}.bam",
        bams = collect_runs_from_sample,
        bams_idx = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_sample(wildcards)))
        #bam_idx = config["results"] + "alignments/merged/{sample}.bam.bai",
    params:
        bams = lambda wildcards, input: list(map(lambda bam: "-I " + bam, input.bams)),
    output:
        vcf = config["results"] + "gvcf/{sample_run}.g.vcf.gz",
        tbi = config["results"] + "gvcf/{sample_run}.g.vcf.gz.tbi",
    conda: "../envs/gatk.yaml"
    threads: 9  # 4 is default
    resources: nodes = 9
    # -I {input.bam} \
    shell: """
        gatk --java-options '-Xmx16g' HaplotypeCaller \
            -R {input.ref} \
            {params.bams} \
            -O {output.vcf} \
            -ERC GVCF \
            --native-pair-hmm-threads {threads} \
            --tmp-dir ~/tmp/{rule}"""

# Consolidation of GVCFs
rule create_sample_map:
    """Create sample map that contains names and paths to all VCFs to be used in consolidate rule.
    Note: This output file will need to be deleted if changing what will be added in the consolidate rule."""
    input:
        gvcfs = expand("{results}gvcf/{sample}.g.vcf.gz",
            results=config["results"],
            sample=SAMPLES),
    output:
        #sample_map = temp(config["results"] + "db/{dataset}.sample_map"),
        sample_map = config["results"] + "db/{dataset}.sample-map",
    threads: 1
    resources: nodes = 1
    run:
        with open(output.sample_map, "w") as sample_map:
            for gvcf in input.gvcfs:
                sample = gvcf.split("/")[-1].split(".")[0]
                sample_map.write(f"{sample}\t{gvcf}\n")

rule list_chromosomes:
    """Find all chromosomes from headers of reference genome.
    
     Each line of the output file is just one chromosome's name.
     The search only keeps numbered chromosomes (not those prefixed with "chr" or any other letters)
     as well as X, Y, and MT. Unplaced contigs are ignored."""
    input:
        ref = config["ref_fasta"],
    output:
        config["results"] + "db/chromosomes.list",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        zcat {input.ref} | grep "^>" | cut -c 2- | cut -d " " -f 1 | grep -x -E "^[0-9]+|X|Y|MT" > {output}
        """

rule consolidate:
    """Combine the chromosomes of .g.vcf files into GenomicsDB datastore."""
    input:
        contigs = config["results"] + "db/chromosomes.list",
        sample_map = config["results"] + "db/{dataset}.sample-map",
    output:
        config["results"] + "db/created_{dataset}.txt",
    params:
        db = config["results"] + "db/{dataset}",
        # Higher value requires more memory and number of file descriptor able to be used at the same time.
        # Attempted with 6, but ran out of memory. Once even 4 was too much. Though did work with 4 when I had fewer samples.
        parallel_intervals = 3,  
    threads: 2  # Just for opening multiple .vcf files at once.
    resources: nodes = 2
    conda: "../envs/gatk.yaml"
    shell: """
        CONTIGS=$(awk 'BEGIN {{ORS = ","}} {{print $0}}' {input.contigs}); \
        if [ -d {params.db} ]; \
        then WORKSPACE_FLAG="genomicsdb-update-workspace-path"; \
        else WORKSPACE_FLAG="genomicsdb-workspace-path"; \
        fi; \
        gatk --java-options '-Xmx8g' GenomicsDBImport \
            --$WORKSPACE_FLAG {params.db} \
            --intervals {input.contigs} \
            --sample-name-map {input.sample_map} \
            --batch-size 50 \
            --genomicsdb-shared-posixfs-optimizations true \
            --reader-threads {threads} \
            --max-num-intervals-to-import-in-parallel {params.parallel_intervals} \
        && touch {output}
        """

## Jointly call variants
rule joint_call_cohort:
    """Use GenomicsDB to jointly call a VCF file."""
    input:
        ref = config["ref_fasta"],
        # Note: The actual text file isn't what is required, but the datastore directory.
        # This .txt file is, however, is created only after the datastore has finished being built.
        db = config["results"] + "db/created_{dataset}.txt",
    output:
        vcf = config["results"] + "joint_call/polyallelic/{dataset}.chr{chr}.vcf.gz",
        #tbi = config["results"] + "joint_call/polyallelic/{dataset}.chr{chr}.vcf.gz.tbi",
    params:
        db = config["results"] + "db/{dataset}",
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    shell: """
        gatk --java-options '-Xmx16g' GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{params.db} \
            -O {output.vcf} \
            -L {wildcards.chr} \
        """

## Split SNPs and indels
rule biallelics_by_mode:
    """Split into SNP- or indel-only .vcf. Then keeps only biallelic sites."""
    input:
        vcf = config["results"] + "joint_call/polyallelic/{dataset}.chr{chr}.vcf.gz",
        tbi = config["results"] + "joint_call/polyallelic/{dataset}.chr{chr}.vcf.gz.tbi",
        ref_fasta = config["ref_fasta"],
    output:
        split = config["results"] + "joint_call/biallelic/{dataset}.{mode}.chr{chr}.vcf.gz",
    params:
        equality = lambda wildcards: "=" if wildcards.mode == "indel" else "!=",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    # 1) Separate multiallelics into different lines
    # 2) Take only SNPs or indels
    # 3) Merge multiallelics back into same lines
    # 4) Keep only biallelics
    # Alternative:
    # bcftools view {input.vcf} \
    # -M2 \
    # -v snps \
    # -Oz \
    # -o {output.split} \
    shell: """
        bcftools norm {input.vcf} \
            -m-any \
            --fasta-ref {input.ref_fasta} \
            -Ou \
        | bcftools view \
            -e'type{params.equality}"snp"' \
            -Ou \
        | bcftools norm \
            -m+any \
            -Ou \
        | bcftools view \
            -M2 \
            -m2 \
            -Oz \
            -o {output.split} \
        """
