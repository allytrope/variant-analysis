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
    output: config["results"] + "ref/ref_genome_original" + ".fna.gz"
    shell: "rsync rsync://" + config["ref"]["address"] + " {output}"

rule bgzip_ref:
    """Convert gzipped ref to bgzipped (a subtype of .gz)."""
    input: config["results"] + "ref/ref_genome_original" + ".fna.gz"
    output: config["results"] + "ref/ref_genome.fna.gz"
    shell: "gunzip -c {input} | bgzip -c > {output}"

rule index_ref:
    """Create BWA index files for reference genome."""
    input: config["results"] + "ref/ref_genome.fna.gz"
    output: multiext(config["results"] + "ref/ref_genome.fna.gz", ".amb", ".ann", ".bwt", ".pac", ".sa")
    shell: "bwa index {input}"

rule index_ref_fai:
    """Create .fai index for reference genome."""
    # Requires bgzipped reference
    input: config["results"] + "ref/ref_genome.fna.gz"
    output: config["results"] + "ref/ref_genome.fna.gz" + ".fai"
    shell: "samtools fqidx {input}"

rule create_ref_dict:
    """Create .dict file for reference genome."""
    input: config["results"] + "ref/ref_genome.fna.gz"
    output: config["results"] + "ref/ref_genome.dict"
    conda: "../envs/gatk.yaml"
    shell: "gatk CreateSequenceDictionary \
            -R {input}"

def sample_path(sample):
    """Find path for corresponding sample name."""
    idx = SAMPLE_NAMES.index(sample)
    return SAMPLE_PATHS[idx] + "/"

# Alignment
rule align:
    """Align .fastq sequence to reference genome."""
    input: ref=config["results"] + "ref/ref_genome.fna.gz",
           idxs=multiext(config["results"] + "ref/ref_genome.fna.gz", ".amb", ".ann", ".bwt", ".pac", ".sa"),
           read_1=lambda wildcards: sample_path("{sample}".format(sample=wildcards.sample)) + "{sample}.R1.fastq.gz",
           read_2=lambda wildcards: sample_path("{sample}".format(sample=wildcards.sample)) + "{sample}.R2.fastq.gz"
    output: config["results"] + "alignments/{sample}.bam"
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
    input: config["results"] + "alignments/{sample}.bam"
    output: config["results"] + "alignments_rg/{sample}.bam"
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
    input: config["results"] + "alignments_rg/{sample}.bam"
    output: config["results"] + "alignments_sorted/{sample}.bam"
    threads: 12
    conda: "../envs/gatk.yaml"
    shell: "gatk SortSam \
                -I {input} \
                -O {output} \
                --SORT_ORDER coordinate \
                --TMP_DIR ~/tmp/{rule}"

# Post-alignment cleanup
# Sometimes gives errors on cluster. But on repeat by themselves (still on cluster), process follows through.
rule mark_duplicates:
    """Mark duplicates with decimal value 1024."""
    input: config["results"] + "alignments_sorted/{sample}.bam"
    output: mdup=config["results"] + "alignments_marked/{sample}.bam",
            metrics=config["results"] + "alignments_marked/{sample}.marked_dup.metrics.txt"
    # Must lower if receive an error about "too many open files". Conversely, can be increased if system allows.
    # To test system capabilites, use `ulimit -a` and look find "open files"
    params: max_open_files=2400
    threads: 12
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' MarkDuplicates \
                -I {input} \
                -O {output.mdup} \
                -M {output.metrics} \
                --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP {params.max_open_files} \
                --TMP_DIR ~/tmp/{rule}"

rule chrom_remapping:
    """Change names of contigs in .vcf according to given mapping.
    For example, can be used to change project contig ids to those of RefSeq.
    Or to change be between form {1, 2, 3,...} and {chr1, chr2, chr3,...}."""
    input: vcf = config["vcf"],
           map=config["map"]
    output: config["results"] + "ref/ref_vcf_remapped.vcf.gz"
    threads: 12
    shell: "bcftools annotate {input.vcf} \
                --rename-chrs {input.map} \
                -o {output} \
                -Oz"

rule index_vcf:
    """Create and index for .vcf file."""
    input: vcf=config["results"] + "ref/ref_vcf_remapped.vcf.gz"
    output: config["results"] + "ref/ref_vcf_remapped.vcf.gz.tbi"
    threads: 12
    conda: "../envs/gatk.yaml"
    shell: "gatk IndexFeatureFile \
            -I {input}"

rule base_recalibration:
    """Create recalibration table for correcting systemic error in base quality scores."""
    input: ref=config["results"] + "ref/ref_genome.fna.gz",
           bam=config["results"] + "alignments_marked/{sample}.bam",
           known_variants=config["results"] + "ref/ref_vcf_remapped.vcf.gz",
           indexed_known_variants=config["results"] + "ref/ref_vcf_remapped.vcf.gz.tbi"
    output: config["results"] + "BQSR/{sample}.BQSR.recal"
    threads: 12
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' BaseRecalibrator \
                -R {input.ref} \
                -I {input.bam} \
                --known-sites {input.known_variants} \
                -O {output}"
                #--tmp-dir ~/tmp/{rule}

rule apply_base_recalibration:
    """Correct systemic error in base quality scores."""
    input: ref=config["results"] + "ref/ref_genome.fna.gz",
           bam=config["results"] + "alignments_marked/{sample}.bam",
           recal=config["results"] + "BQSR/{sample}.BQSR.recal"
    output: config["results"] + "alignments_recalibrated/{sample}.bam"
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' ApplyBQSR \
                -R {input.ref} \
                -I {input.bam} \
                --bqsr-recal-file {input.recal} \
                -O {output}"

# Variant calling
rule call_variants:
    """Call variants to make VCF file."""
    input: ref=config["results"] + "ref/ref_genome.fna.gz",
           bam=config["results"] + "alignments_recalibrated/{sample}.bam"
    output: config["results"] + "vcf/{sample}.g.vcf.gz"
    threads: 12
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' HaplotypeCaller \
                -R {input.ref} \
                -I {input.bam} \
                -O {output} \
                -ERC GVCF"

rule create_sample_map:
    """Create sample map that contains names and paths to all VCFs to be used in consolidate rule."""
    input: expand("{results}vcf/{sample}.g.vcf.gz", results=config["results"], sample=SAMPLE_NAMES)
    output: config["results"] + "/cohort.sample_map"
    run:
        import os
        with open("cohort.sample_map", "w") as out:
            for path, sample in FULL_SAMPLE_PATHS:
                string = sample + "\t" + config["results"] + "vcf/" + sample + ".g.vcf.gz" + "\n"
                out.write(string)

# Consolidation of GVCFs
rule consolidate:
    """Combines GVCFs into GenomicsDB (a datastore)."""
    input: expand("{results}{sample}.g.vcf", results=config["results"], sample=SAMPLE_NAMES), config["results"] + "cohort.sample_map"
    output: config["results"] + "db/genomicsdb"
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' GenomicsDBImport \
                --sample-name-map cohort.sample_map \
                --genomicsdb-workspace-path db \
                --genomicsdb-update-workspace-path db"

rule joint_call_cohort:
    """Use GenomicsDB to jointly call a VCF file."""
    input: ref=config["results"] + "ref/ref_genome.fna.gz",
           db=config["results"] + "db"
    output: config["results"] + "joint-call/joint-call.vcf.gz"
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' GenotypeGVCFs \
                -R {input.ref} \
                -V gendb://{input.db} \
                -O {output}"

# Variant quality score recalibration
rule calls_recalibration:
    """Build a recalibration model."""
    input: ref=config["results"] + "ref/ref_genome.fna.gz",
           vcf=config["results"] + "joint-call/joint-call.vcf.gz",
           truth=config["VQSR"]["truth_vcf"],
           training=config["VQSR"]["training_vcf"]
    output: recal=config["results"] + "VQSR/VQSR.recal",
            tranches=config["results"] + "VQSR/output.tranches"
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
    input: ref=config["results"] + "ref/ref_genome.fna.gz",
           vcf=config["results"] + "joint-call/joint-call.vcf.gz",
           recal=config["results"] + "VQSR/VQSR.recal",
           tranches=config["results"] + "VQSR/output.tranches"
    output: config["results"] + "joint-call/recalibrated_joint-call.vcf.gz"
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options 'Xmx8g' ApplyVQSR \
                -R {input.ref} \
                -V {input.vcf} \
                -O {output} \
                --recal-file {input.recal} \
                --tranches-file {input.tranches}"

# CHANGE TO WORK FOR MMUL_10 DATA
rule snp_summary:
    """Gives a summary of SNP data from .vcf file."""
    input: vcf=config["results"] + "joint-call/recalibrated_joint-call.vcf.gz"
    output: config["results"] + "joint-call/snp_summary.txt"
    shell: "vep \
    --custom /home/flow/ensembl-vep/Mmul_8.0.1.92.chr.gff3.gz,,gff \
    --fasta Mmul_8.0.1.chromosome.fa \
    --gff /home/flow/ensembl-vep/Mmul_8.0.1.92.chr.gff3.gz \
    --input_file {input.vcf} \
    --output_file {output} \
    --stats_text "