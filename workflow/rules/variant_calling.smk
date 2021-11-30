"""Contain workflow for performing variant analysis from paired-end FASTQ files."""

import gzip
import os

# Path and sample names for input .fastq files
if config["data"]:
    FULL_SAMPLE_PATHS = [os.path.split(path) for path in config["data"]]
    SAMPLE_PATHS, SAMPLE_NAMES = zip(*FULL_SAMPLE_PATHS)
else:
    FULL_SAMPLE_PATHS = SAMPLE_PATHS = SAMPLE_NAMES = []

# Prepare reference genome
rule download_ref:
    """Download reference genome."""
    output: config["ref"]["fasta"]
    shell: "rsync rsync://" + config["ref"]["address"] + " {output}"

rule index_ref:
    """Create BWA index files for reference genome."""
    input: config["ref"]["fasta"]
    output: multiext(config["ref"]["fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    shell: "bwa index {input}"

def sample_path(sample):
    """Find path for corresponding sample name."""
    idx = SAMPLE_NAMES.index(sample)
    return SAMPLE_PATHS[idx] + "/"

# Alignment
rule align:
    """Align .fastq sequence to reference genome."""
    input: ref=config["ref"]["fasta"],
           idxs=multiext(config["ref"]["fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
           read_1=lambda wildcards: sample_path("{sample}".format(sample=wildcards.sample)) + "{sample}.R1.fastq.gz",
           read_2=lambda wildcards: sample_path("{sample}".format(sample=wildcards.sample)) + "{sample}.R2.fastq.gz"
    output: config["output_path"] + "alignments/{sample}.bam"
    threads: 24
    shell: "bwa mem -t {threads} {input.ref} {input.read_1} {input.read_2} | \
            samtools view -Sb - > {output}"

def find_read_group_info(sample):
    """Find corresponding .fastq file from sample name. For Casava 1.8 Illumina header implementation."""
    for file in config["data"]:
        if file.endswith(sample):
            with gzip.open(file + ".R1.fastq.gz", "rt") as fastq:
                line = fastq.readline()
                header = line.split(":")
                rgid = header[2] + "." + header[3]
                rgpu = header[2] + "." + header[3] + "." + sample
                rgsm = sample
                rgpl = "ILLUMINA"
                rglb = sample  # Using sample name to make every .fastq file count as separate library
                print("Read group info:", rgid, rgpu, rgsm, rgpl, rglb)
                return rgid, rgpu, rgsm, rgpl, rglb

rule add_read_groups:
    """Add @RG row to .bam header."""
    input: config["output_path"] + "alignments/{sample}.bam"
    output: config["output_path"] + "alignments_rg/{sample}.bam"
    params: header_rg = lambda wildcards: find_read_group_info("{sample}".format(sample=wildcards.sample))
    conda: "../envs/gatk.yaml"
    shell: "gatk AddOrReplaceReadGroups \
      -I {input} \
      -O {output} \
      -RGID {params.header_rg[0]} \
      -RGPU {params.header_rg[1]} \
      -RGSM {params.header_rg[2]} \
      -RGPL {params.header_rg[3]} \
      -RGLB {params.header_rg[4]} \
      "

rule sort_reads:
    """Sort reads in .bam file."""
    input: config["output_path"] + "alignments_rg/{sample}.bam"
    output: config["output_path"] + "sorted_alignments/{sample}.bam"
    conda: "../envs/gatk.yaml"
    shell: "gatk SortSam \
                -I {input} \
                -O {output} \
                -SORT_ORDER coordinate"

# Post-alignment cleanup
rule mark_duplicates:
    """Mark duplicates with decimal value 1024."""
    input: config["output_path"] + "sorted_alignments/{sample}.bam"
    output: mdup=config["output_path"] + "marked_bam/{sample}.bam",
            metrics=config["output_path"] + "marked_bam/{sample}.marked_dup.metrics.txt"
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' MarkDuplicates \
                -I {input} \
                -O {output.mdup} \
                -M {output.metrics}"

rule base_recalibration:
    """Create recalibration table for correcting systemic error in base quality scores."""
    input: ref=config["ref"]["fasta"],
           bam=config["output_path"] + "marked_bam/{sample}.bam",
           known_variants=config["vcf"]
    output: config["output_path"] + "BQSR/{sample}.BQSR.recal"
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' BaseRecalibrator \
                -R {input.ref} \
                -I {input.bam} \
                --known-sites {input.known_variants} \
                -O {output}"

rule apply_base_recalibration:
    """Correct systemic error in base quality scores."""
    input: ref=config["ref"]["fasta"],
           bam=config["output_path"] + "marked-bam/{sample}.bam",
           recal=config["output_path"] + "BQSR/{sample}.BQSR.recal"
    output: config["output_path"] + "cleaned_bam/{sample}.bam"
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' ApplyBQSR \
                -R {input.ref} \
                -I {input.bam} \
                --bqsr-recal-file {input.recal} \
                -O {output}"

# Variant calling
rule call_variants:
    """Call variants to make VCF file."""
    input: ref=config["ref"]["fasta"],
           bam=config["output_path"] + "cleaned_bam/{sample}.bam"
    output: config["output_path"] + "vcf/{sample}.g.vcf.gz"
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' HaplotypeCaller \
                -R {input.ref} \
                -I {input.bam} \
                -O {output} \
                -ERC GVCF"

rule create_sample_map:
    """Create sample map that contains names and paths to all VCFs to be used in consolidate rule."""
    input: expand("{output_path}vcf/{sample}.g.vcf.gz", output_path=config["output_path"], sample=SAMPLE_NAMES)
    output: config["output_path"] + "/cohort.sample_map"
    run:
        import os
        with open("cohort.sample_map", "w") as out:
            for path, sample in FULL_SAMPLE_PATHS:
                string = sample + "\t" + config["output_path"] + "vcf/" + sample + ".g.vcf.gz" + "\n"
                out.write(string)

# Consolidation of GVCFs
rule consolidate:
    """Combines GVCFs into GenomicsDB (a datastore)."""
    input: expand("{output_path}{sample}.g.vcf", output_path=config["output_path"], sample=SAMPLE_NAMES), config["output_path"] + "cohort.sample_map"
    output: config["output_path"] + "db/genomicsdb"
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' GenomicsDBImport \
                --sample-name-map cohort.sample_map \
                --genomicsdb-workspace-path db \
                --genomicsdb-update-workspace-path db"

rule joint_call_cohort:
    """Use GenomicsDB to jointly call a VCF file."""
    input: ref=config["ref"]["fasta"],
           db=config["output_path"] + "db"
    output: config["output_path"] + "joint-call/joint-call.vcf.gz"
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' GenotypeGVCFs \
                -R {input.ref} \
                -V gendb://{input.db} \
                -O {output}"

# Variant quality score recalibration
rule calls_recalibration:
    """Build a recalibration model."""
    input: ref=config["ref"]["fasta"],
           vcf=config["output_path"] + "joint-call/joint-call.vcf.gz",
           truth=config["VQSR"]["truth_vcf"],
           training=config["VQSR"]["training_vcf"]
    output: recal=config["output_path"] + "VQSR/VQSR.recal",
            tranches=config["output_path"] + "VQSR/output.tranches"
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' VariantRecalibrator \
                -R {input.ref} \
                -V {input.vcf} \
                -O {output.recal} \
                --resource source1,known=false,truth=true,training=true,prior=10.0:{input.truth} \
                --resource source2,known=false,truth=false,training=true,prior=10.0:{input.training} \
                #-an add anotations from VCF \
                --tranches-file {output.tranches}"

rule apply_variant_recalibration:
    """Apply recalibration model."""
    input: ref=config["ref"]["fasta"],
           vcf=config["output_path"] + "joint-call/joint-call.vcf.gz",
           recal=config["output_path"] + "VQSR/VQSR.recal",
           tranches=config["output_path"] + "VQSR/output.tranches"
    output: config["output_path"] + "joint-call/recalibrated_joint-call.vcf.gz"
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options 'Xmx8g' ApplyVQSR \
                -R {input.ref} \
                -V {input.vcf} \
                -O {output} \
                --recal-file {input.recal} \
                --tranches-file {input.tranches}"
                