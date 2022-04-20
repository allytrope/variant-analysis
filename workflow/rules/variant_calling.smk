"""Contain workflow for performing variant analysis from paired-end FASTQ files."""

import gzip
import os

# Path and sample names for input .fastq files
if config["data"]:
    FULL_SAMPLE_PATHS = [os.path.split(path) for path in config["data"]]
    SAMPLE_PATHS, SAMPLE_NAMES = zip(*FULL_SAMPLE_PATHS)
else:
    FULL_SAMPLE_PATHS = SAMPLE_PATHS = SAMPLE_NAMES = []

CONFIG = config["variant_calling"]

# Prepare reference genome
rule download_ref:
    """Download reference genome."""
    output: config["results"] + "ref/ref_genome_original" + ".fna.gz",
    shell: "rsync rsync://" + config["ref"]["address"] + " {output}"

rule bgzip_ref:
    """Convert gzipped ref to bgzipped (a subtype of .gz)."""
    input: config["results"] + "ref/ref_genome_original" + ".fna.gz",
    output: config["results"] + "ref/ref_genome.fna.gz",
    shell: "gunzip -c {input} | bgzip -c > {output}"

rule index_ref:
    """Create BWA index files for reference genome."""
    input: config["results"] + "ref/ref_genome.fna.gz",
    output: multiext(config["results"] + "ref/ref_genome.fna.gz", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    conda: "../envs/bio.yaml"
    shell: "bwa index {input}"

rule index_ref_fai:
    """Create .fai index for reference genome."""
    # Requires bgzipped reference
    input: config["results"] + "ref/ref_genome.fna.gz",
    output: config["results"] + "ref/ref_genome.fna.gz.fai",
    conda: "../envs/bio.yaml"
    shell: "samtools fqidx {input}"

rule create_ref_dict:
    """Create .dict file for reference genome."""
    input: config["results"] + "ref/ref_genome.fna.gz",
    output: config["results"] + "ref/ref_genome.dict",
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
    input: ref = config["results"] + "ref/ref_genome.fna.gz",
           idxs = multiext(config["results"] + "ref/ref_genome.fna.gz", ".amb", ".ann", ".bwt", ".pac", ".sa"),
           read_1 = lambda wildcards: sample_path("{sample}".format(sample=wildcards.sample)) + "{sample}.R1.fastq.gz",
           read_2 = lambda wildcards: sample_path("{sample}".format(sample=wildcards.sample)) + "{sample}.R2.fastq.gz",
    output: config["results"] + "alignments/{sample}.bam",
    threads: 24
    conda: "../envs/bio.yaml"
    shell: "bwa mem -t {threads} {input.ref} {input.read_1} {input.read_2} | \
            samtools view -Sb - > {output}"

'''
rule align_with_rg:
    """Align .fastq sequence to reference genome."""
    input: ref=config["results"] + "ref/ref_genome.fna.gz",
           idxs=multiext(config["results"] + "ref/ref_genome.fna.gz", ".amb", ".ann", ".bwt", ".pac", ".sa"),
           read_1=lambda wildcards: sample_path("{sample}".format(sample=wildcards.sample)) + "{sample}.R1.fastq.gz",
           read_2=lambda wildcards: sample_path("{sample}".format(sample=wildcards.sample)) + "{sample}.R2.fastq.gz"
    output: config["results"] + "alignments/{sample}.bam"
    threads: 24
    conda: "../envs/bio.yaml"
    shell: "bwa mem {input.ref} {input.read_1} {input.read_2} \
                -t {threads} \
                -R '' | \
            samtools view -Sb - > {output}"
'''

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
    input: config["results"] + "alignments/{sample}.bam",
    output: config["results"] + "alignments_rg/{sample}.bam",
    params: header_rg = lambda wildcards: find_read_group_info("{sample}".format(sample=wildcards.sample)),
    threads: 1
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
    input: config["results"] + "alignments_rg/{sample}.bam",
    output: config["results"] + "alignments_sorted/{sample}.bam",
    threads: 1
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
    input: config["results"] + "alignments_sorted/{sample}.bam",
    output: mdup = config["results"] + "alignments_marked/{sample}.bam",
            metrics = config["results"] + "alignments_marked/{sample}.marked_dup.metrics.txt",
    # Must lower if receive an error about "too many open files". Conversely, can be increased if system allows.
    # To test system capabilites, use `ulimit -a` and look find "open files"
    params: max_open_files = 2400,
    threads: 1
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' MarkDuplicates \
                -I {input} \
                -O {output.mdup} \
                -M {output.metrics} \
                --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP {params.max_open_files} \
                --TMP_DIR ~/tmp/{rule}"

rule contig_remapping:
    """Change names of contigs in .vcf according to given mapping.
    For example, can be used to change project contig ids to those of RefSeq.
    Or to change be between form {1, 2, 3,...} and {chr1, chr2, chr3,...}."""
    input: vcf = CONFIG["vcf"],
           map = CONFIG["contig_remapping"],
    output: config["results"] + "ref/ref_vcf_remapped.vcf.gz",
    threads: 1
    conda: "../envs/bio.yaml"
    shell: "bcftools annotate {input.vcf} \
                --rename-chrs {input.map} \
                -o {output} \
                -Oz"

rule tbi_index:
    """Create .tbi index."""
    input: "{path}/{name}.vcf.gz",
    output: "{path}/{name}.vcf.gz.tbi",
    threads: 1
    conda: "../envs/gatk.yaml"
    shell: "gatk IndexFeatureFile \
            -I {input}"

rule base_recalibration:
    """Create recalibration table for correcting systemic error in base quality scores."""
    input: ref = config["results"] + "ref/ref_genome.fna.gz",
           bam = config["results"] + "alignments_marked/{sample}.bam",
           known_variants = CONFIG["BQSR_known_variants"],
           indexed_known_variants = CONFIG["BQSR_known_variants"] + ".tbi",
    output: config["results"] + "BQSR/{sample}.BQSR.recal",
    threads: 1
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' BaseRecalibrator \
                -R {input.ref} \
                -I {input.bam} \
                --known-sites {input.known_variants} \
                -O {output}"
                #--tmp-dir ~/tmp/{rule}

rule apply_base_recalibration:
    """Correct systemic error in base quality scores."""
    input: ref = config["results"] + "ref/ref_genome.fna.gz",
           bam = config["results"] + "alignments_marked/{sample}.bam",
           recal = config["results"] + "BQSR/{sample}.BQSR.recal",
    output: config["results"] + "alignments_recalibrated/{sample}.bam",
    threads: 1
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' ApplyBQSR \
                -R {input.ref} \
                -I {input.bam} \
                --bqsr-recal-file {input.recal} \
                -O {output}"

# Variant calling
rule call_variants:
    """Call variants to make VCF file."""
    input: ref = config["results"] + "ref/ref_genome.fna.gz",
           bam = config["results"] + "alignments_recalibrated/{sample}.bam",
    output: vcf = config["results"] + "vcf/{sample}.g.vcf.gz",
            tbi = config["results"] + "vcf/{sample}.g.vcf.gz.tbi",  # Index file
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
    input: gvcfs = expand("{results}vcf/{sample}.g.vcf.gz", results=config["results"], sample=SAMPLE_NAMES),
    output: sample_map = config["results"] + "ref/{workspace}.sample_map",
    run:
        import os
        with open(output.sample_map, "w") as out:
            for path, sample in FULL_SAMPLE_PATHS:
                string = sample + "\t" + config["results"] + "vcf/" + sample + ".g.vcf.gz" + "\n"
                out.write(string)

rule find_placed_contigs:
    """Find all contigs from reference that start with 'NC_'.
    This is the RefSeq designation for placed contigs (chromosomes + mitochodrion).
    If chromosomes use a different naming system, change or remove the `grep` command."""
    input: config["results"] + "ref/ref_genome.fna.gz",
    output: config["results"] + "ref/placed_contigs.list",
    threads: 1
    conda: "../envs/bio.yaml"
    shell: "samtools view {input} | cut -f 1 | grep '^NC_' > {output}"

rule consolidate:
    """Combine the placed_contigs of .g.vcf files into GenomicsDB datastore."""
    input: contigs = config["results"] + "ref/placed_contigs.list",
           sample_map = config["results"] + "ref/{workspace}.sample_map",
           #gvcfs=expand("{results}vcf/{sample}.g.vcf.gz", results=config["results"], sample=SAMPLE_NAMES)
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

rule joint_call_cohort:
    """Use GenomicsDB to jointly call a VCF file."""
    input: ref = config["results"] + "ref/ref_genome.fna.gz",
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
    input: ref = config["results"] + "ref/ref_genome.fna.gz",
           db = config["results"] + "db/{workspace}",
    output: config["results"] + "joint_call/{workspace}.{chrom}.vcf.gz",
    threads: 1
    conda: "../envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx8g' GenotypeGVCFs \
                -R {input.ref} \
                -V gendb://{input.db} \
                -O {output} \
                -L {wildcards.chrom}"
                
# Copied from above rule
# rule chromosome_remapping_per_chromosome:
#     """Change names of contigs in .vcf according to given mapping.
#     For example, can be used to change project contig ids to those of RefSeq.
#     Or to change be between form {1, 2, 3,...} and {chr1, chr2, chr3,...}."""
#     input: vcf = "{path}/{name}."
#            map = CONFIG["contig_remapping"],
#     output: config["results"] + "ref/ref_vcf_remapped.vcf.gz",
#     threads: 1
#     conda: "../envs/bio.yaml"
#     shell: "bcftools annotate {input.vcf} \
#                 --rename-chrs {input.map} \
#                 -o {output} \
#                 -Oz"