"""Contain workflow for performing variant analysis from paired-end FASTQ files."""

import gzip
import os


CONFIG = config["variant_calling"]

## Creating indicies
rule index_ref:
    """Create BWA index files for reference genome."""
    input: CONFIG["ref_fasta"],
    output: multiext(CONFIG["ref_fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    conda: "../envs/bio.yaml"
    shell: "bwa index {input}"

#For some reason is causing a JSONDecodeError
# rule index_ref_fai:
#     """Create .fai index for reference genome."""
#     # Requires bgzipped reference
#     input: CONFIG["ref_fasta"],
#     output: CONFIG["ref_fasta"] + ".fai",
#     conda: "../envs/bio.yaml"
#     shell: "samtools fqidx {input}"

rule create_ref_dict:
    """Create .dict file for reference genome."""
    #wildcard_constraints: fasta = "fasta|fna|fa"
    input: fasta = CONFIG["ref_fasta"],
    output: dict = ".".join(CONFIG["ref_fasta"].split(".")[:-2]) + ".dict",  # Replaces the ".fna.gz" ending with ".dict"
    conda: "../envs/gatk.yaml"
    shell: "gatk CreateSequenceDictionary \
            -R {input.fasta}"

rule tbi_index:
    """Create .tbi index."""
    input: "{path}/{name}.vcf.gz",
    output: "{path}/{name}.vcf.gz.tbi",
    threads: 1
    conda: "../envs/bio.yaml"
    shell: "bcftools index {input} \
                --tbi"

rule bai_index:
    """Create .bai index for .bam file."""
    input: "{path}/{sample}.bam",
    output: "{path}/{sample}.bam.bai",
    conda: "../envs/bio.yaml"
    shell: "samtools index {input}"

## Trimming
rule cut_adaptors:
    """Cut out 3' end adapters."""
    input: #reads = CONFIG["reads"]
           read1 = config["resources"] + "reads/{sample}.R1.fastq.gz",
           read2 = config["resources"] + "reads/{sample}.R2.fastq.gz",
           adapters = CONFIG["adapters"],
    output: trimmed1 = config["results"] + "trimmed/{sample}.trimmed.R1.fastq.gz",
            trimmed2 = config["results"] + "trimmed/{sample}.trimmed.R2.fastq.gz",
    threads: 4
    conda: "../envs/bio.yaml"
    shell: """
        R1_ADAPTER=$(awk '$1 == "{wildcards.sample}" {{print $2}}' {input.adapters}); \
        R2_ADAPTER=$(awk '$1 == "{wildcards.sample}" {{print $3}}' {input.adapters}); \
        cutadapt {input.read1} {input.read2} \
            -a $R1_ADAPTER \
            -A $R2_ADAPTER \
            --cores {threads} \
            -o {output.trimmed1} \
            -p {output.trimmed2}"""

## Alignment
rule align_with_rg:
    """Align .fastq sequence to reference genome.

    Extracts read group information. Then aligns to references genome while adding this read group info.
    Then sorts the output by coordinates."""
    input: ref = CONFIG["ref_fasta"],
           idxs = multiext(CONFIG["ref_fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
           trimmed1 = config["results"] + "trimmed/{sample}.trimmed.R1.fastq.gz",
           trimmed2 = config["results"] + "trimmed/{sample}.trimmed.R2.fastq.gz",
    output: alignment = config["results"] + "alignments_raw/{sample}.bam",
    threads: 24
    conda: "../envs/bio.yaml"
    # First of RG's tags must be SM and last must be PU because of how I have to call the sample names.
    shell: """
        bash workflow/scripts/align.sh {input.trimmed1} {input.trimmed2} {input.ref} {output} {threads} {wildcards.sample}"""

rule alignment_QC:
    input: alignment = config["results"] + "alignments_raw/{sample}.bam",
    output: alignment = config["results"] + "alignments/{sample}.bam",
    conda: "../envs/bio.yaml"
    shell: """
        samtools sort {input.alignment} \
            -n \
        | samtools fixmate - - \
            -m \
        | samtools sort - \
        | samtools markdup - {output.alignment} \
            -T ~/tmp/{rule}
        """

## Base recalibration
rule base_recalibration:
    """Create recalibration table for correcting systemic error in base quality scores."""
    input: ref = CONFIG["ref_fasta"],
           ref_idx = CONFIG["ref_fasta"] + ".fai",
           ref_dict = ".".join(CONFIG["ref_fasta"].split(".")[:-2]) + ".dict",
           bam = config["results"] + "alignments/{sample}.bam",
           known_variants = CONFIG["BQSR_known_variants"],
           indexed_known_variants = CONFIG["BQSR_known_variants"] + ".tbi",
    output: config["results"] + "BQSR/{sample}.BQSR.recal",
    threads: 1
    conda: "../envs/gatk.yaml"
    # Note tmp directory must already exist
    shell: "gatk --java-options '-Xmx8g' BaseRecalibrator \
                -R {input.ref} \
                -I {input.bam} \
                --known-sites {input.known_variants} \
                -O {output} \
                --tmp-dir ~/tmp/{rule}"

rule apply_base_recalibration:
    """Correct systemic error in base quality scores."""
    input: ref = CONFIG["ref_fasta"],
           ref_idx = CONFIG["ref_fasta"] + ".fai",
           ref_dict = ".".join(CONFIG["ref_fasta"].split(".")[:-2]) + ".dict",
           bam = config["results"] + "alignments/{sample}.bam",
           recal = config["results"] + "BQSR/{sample}.BQSR.recal",
    output: config["results"] + "alignments_recalibrated/{sample}.bam",
    threads: 1
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' ApplyBQSR \
                -R {input.ref} \
                -I {input.bam} \
                --bqsr-recal-file {input.recal} \
                -O {output}"

## Variant calling
# Modified input to skip recalibration
rule call_variants:
    """Call variants to make VCF file."""
    input: ref = CONFIG["ref_fasta"],
           bam = config["results"] + "alignments/{sample}.bam",
           bam_idx = config["results"] + "alignments/{sample}.bam.bai",
    output: vcf = config["results"] + "gvcf/{sample}.g.vcf.gz",
            tbi = config["results"] + "gvcf/{sample}.g.vcf.gz.tbi",  # Index file
    conda: "../envs/gatk.yaml"
    threads: 24  # 4 is default
    shell: "gatk --java-options '-Xmx16g' HaplotypeCaller \
                -R {input.ref} \
                -I {input.bam} \
                -O {output.vcf} \
                -ERC GVCF \
                --native-pair-hmm-threads {threads}"

# Consolidation of GVCFs
rule create_sample_map:
    """Create sample map that contains names and paths to all VCFs to be used in consolidate rule.
    Note: This output file will need to be deleted if changing what will be added in the consolidate rule."""
    input: gvcfs = expand("{results}gvcf/{sample}.g.vcf.gz", results=config["results"], sample=SAMPLE_NAMES),
    output: sample_map = config["results"] + "db/{workspace}.sample_map",
    params: in_path = config["results"] + "gvcf/",
    threads: 1
    shell: """
        ls {params.in_path} | awk -v FS='.' -v OFS='\t' '/gz$/ {{print $1,"{params.in_path}"$0}}' > {output.sample_map}
        """

rule find_chromosomes:
    """Find all chromosomes from reference genome. Each line of the output file is just one chromosome's name."""
    input: CONFIG["ref_fasta"],  # Or something different
    output: config["results"] + "db/chromosomes.list",
    threads: 1
    conda: "../envs/bio.yaml"
    # Need to generalize for whatever the count of chromosomes are.
    shell: """
        echo {{1..22}} | awk -v RS=' ' '{{print $1}} END {{print "X\\nY\\nMT"}}' > {output}
        """

rule consolidate:
    """Combine the chromosomes of .g.vcf files into GenomicsDB datastore."""
    input: contigs = config["results"] + "db/chromosomes.list",
           sample_map = config["results"] + "db/{workspace}.sample_map",
    output: config["results"] + "db/created_{workspace}.txt",
    params: db = config["results"] + "db/{workspace}",
            parallel_intervals = 4,  # Higher value requires more memory and number of file descriptor able to be used at the same time
    threads: 2  # Just for opening multiple .vcf files at once.
    conda: "../envs/gatk.yaml"
    shell:  """
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
                --max-num-intervals-to-import-in-parallel {params.parallel_intervals}; \
            touch {output}
            """

## Jointly call variants
rule joint_call_cohort:
    """Use GenomicsDB to jointly call a VCF file."""
    input: ref = CONFIG["ref_fasta"],
           db = config["results"] + "db/{workspace}",
    output: config["results"] + "joint_call/{workspace}.vcf.gz",
    threads: 1
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' GenotypeGVCFs \
                -R {input.ref} \
                -V gendb://{input.db} \
                -O {output}"

rule joint_call_cohort_per_chromosome:
    """Use GenomicsDB to jointly call a VCF file."""
    input: ref = CONFIG["ref_fasta"],
           db = config["results"] + "db/{workspace}",
    output: config["results"] + "joint_call/{workspace}.{chrom}.vcf.gz",
    threads: 1
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' GenotypeGVCFs \
                -R {input.ref} \
                -V gendb://{input.db} \
                -O {output} \
                -L {wildcards.chrom}"

## Split SNPs and indels
rule subset_mode:
    """Split into SNP- or indel-only .vcf. The wildcard `mode` can be "SNP" or "indel". Then keeps only biallelic sites."""
    input: vcf = config["results"] + "joint_call/{workspace}.vcf.gz",
           vcf_index = config["results"] + "joint_call/{workspace}.vcf.gz.tbi",
           ref_fasta = config["variant_calling"]["ref_fasta"],
    output: split = config["results"] + "joint_call/split/{workspace}.{mode}.biallelic.vcf.gz",
    params: equality = lambda wildcards: "=" if wildcards.mode == "indel" else "!=",
    threads: 1
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
            -Oz \
            -o {output.split} \
        """

