"""Contain workflow for performing variant analysis from paired-end FASTQ files."""

from Bio.Seq import Seq
import gzip
import os
import pandas as pd


## Trimming
rule cut_adapters_with_i5_i7:
    """Cut out 3' end adapters. Uses i7 and i5 adapter information. This implementation takes information from the first line of the .R1.fastq file.
    This thus assumes that all reads in both files use the same i7 and i5 adapters. Additionally, this assumes that the Truseq Dual Index Library was used for the adapters
    surrounding the i7 and i5 adapters. These are hard coded under "params". The i7 information is used for read 1 and i5 for read 2.
    Additionally, if the i7 and i5 adapters must actually be in the file headers."""
    wildcard_constraints:
        seq = "WGS|WES"
    input: 
        reads = lambda wildcards: expand(config["reads"] + "{seq}{sample}.{read}.fastq.gz",
            seq=wildcards.seq,
            sample=wildcards.sample,
            read=["R1", "R2"]),
    output: 
        trimmed = expand(config["results"] + "trimmed/{{seq}}{{sample}}.{read}.fastq.gz",
            read=["R1", "R2"]),
    params: 
        pre_i7 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        post_i7 = "ATCTCGTATGCCGTCTTCTGCTTG",
        pre_i5 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        post_i5 = "GTGTAGATCTCGGTGGTCGCCGTATCATT",
    threads: 4
    resources: nodes = 4
    conda: "../envs/bio.yaml"
    shell: """
        echo $SHELL;
        ADAPTERS=$(gunzip -c {input.reads[0]} | head -n 1 | cut -d " " -f 2 | cut -d ":" -f 4);
        i7=$(echo $ADAPTERS | cut -d "+" -f 1);
        i5=$(echo $ADAPTERS | cut -d "+" -f 2);
        R1_END_ADAPTER="{params.pre_i7}${{i7}}{params.post_i7}";
        R2_END_ADAPTER="{params.pre_i5}${{i5}}{params.post_i5}";
        cutadapt {input.reads} \
            -a $R1_END_ADAPTER \
            -A $R2_END_ADAPTER \
            --cores {threads} \
            -o {output.trimmed[0]} \
            -p {output.trimmed[1]}"""

# Unfinished
def find_barcode(wildcards, input):
    """Find sample barcode from barcodes file. Barcode file contains reverse complements."""
    df = pd.read_table(input.barcodes, header=None, names=["sample", "barcodes"])
    barcode = df.set_index("sample")["barcodes"].to_dict()[wildcards.sample]
    return str(Seq(barcode).reverse_complement)
rule cut_GBS_adapters:
    """Cut out 3' end adapters. This way differs from the other adapter removing rule.

    Used for the GBS that I have. With this data, there is no i7 or i5 values given in the header lines.
    Instead, R1 and R2 reads have the adapters given under params.
    After both of these adapters, there is a pseudoconsensus sequence, but left out since Cutadapt will remove anything after the given adapters.
    However, before the R2 is a per sample consenus sequence that must be discovered"""
    wildcard_constraints:
        seq = "GBS",
    input:
        reads = lambda wildcards: expand(config["reads"] + "{seq}{sample}.{read}.fastq.gz",
            seq=wildcards.seq,
            sample=wildcards.sample,
            read=["R1", "R2"]),
        barcodes = config["barcodes"],
    output: 
        trimmed = expand(config["results"] + "trimmed/{{seq}}{{sample}}.{read}.fastq.gz",
            read=["R1", "R2"]),
    params: 
        barcode = lambda wildcards, input: find_barcode(wildcards, input),
        R1_adapter = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",
        R2_adapter = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    threads: 4
    resources: nodes = 4
    conda: "../envs/bio.yaml"
    shell: """
        echo $SHELL;
        ADAPTERS=$(gunzip -c {input.reads[0]} | head -n 1 | cut -d " " -f 2 | cut -d ":" -f 4);
        R1_END_ADAPTER="{params.R1_adapter}";
        R2_END_ADAPTER="{params.barcode}{params.R2_adapter}";
        cutadapt {input.reads} \
            -a $R1_END_ADAPTER \
            -A $R2_END_ADAPTER \
            --cores {threads} \
            -o {output.trimmed[0]} \
            -p {output.trimmed[1]}"""

## Alignment
rule align:
    """Align .fastq sequence to reference genome.

    Extracts read group information. Then aligns to references genome while adding this read group info.
    Then sorts the output by coordinates."""
    input:
        ref = config["ref_fasta"],
        idxs = multiext(config["ref_fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
        trimmed = lambda wildcards: expand(config["results"] + "trimmed/{seq}{sample}.{read}.fastq.gz",
            seq=wildcards.seq,
            sample=wildcards.sample,
            read=["R1", "R2"]),
    output:
        alignment = config["results"] + "alignments/raw/{seq}{sample}.bam",
    threads: 24
    resources: nodes = 24
    conda: "../envs/bio.yaml"
    # First of RG's tags must be SM and last must be PU because of how I have to call the sample names.
    shell: """
        bash workflow/scripts/align.sh {input.trimmed[0]} {input.trimmed[1]} {input.ref} {output} {threads} {wildcards.sample}"""

rule alignment_markdup:
    """Fix mate pairs and mark duplicate reads."""
    input:
        alignment = config["results"] + "alignments/raw/{seq}{sample}.bam",
    output:
        alignment = config["results"] + "alignments/markdup/{seq}{sample}.bam",
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
            -T ~/tmp/{rule}
        """

# Not currently used.
# rule merge_like_samples_plus_markdup:
#     """Combine .bam files from multiple runs of same sample."""
#     input:
#         split_alignments = lambda wildcards: expand(config["results"] + "alignments/markdup/{seq}{sample}.bam",
#             sample=wildcards.sample,
#             run=RUNS[wildcards.sample]),
#     output:
#         merged_alignment = config["results"] + "alignments/merged/{sample}.bam",
#     threads: 1
#     resources: nodes = 1
#     shell: """
#         samtools merge {input.split_alignments} \
#             -o {output.merged_alignment}
#         """

## Base recalibration
rule base_recalibration:
    """Create recalibration table for correcting systemic error in base quality scores."""
    input: 
        ref = config["ref_fasta"],
        ref_idx = config["ref_fasta"] + ".fai",
        ref_dict = ".".join(config["ref_fasta"].split(".")[:-2]) + ".dict",
        bam = config["results"] + "alignments/markdup/{sample}.bam",
        known_variants = config["BQSR_known_variants"],
        indexed_known_variants = config["BQSR_known_variants"] + ".tbi",
    output:
        config["results"] + "alignments/recalibrated/recal_tables/{sample}.BQSR.recal",
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
        bam = config["results"] + "alignments/{sample}.bam",
        recal = config["results"] + "alignments/recalibrated/recal_tables/{sample}.BQSR.recal",
    output:
        config["results"] + "alignments/recalibrated/{sample}.bam",
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    # Telling not to make bam index since it automatically writes as 
    shell: """
        gatk --java-options '-Xmx8g' ApplyBQSR \
            -R {input.ref} \
            -I {input.bam} \
            --bqsr-recal-file {input.recal} \
            -O {output} \
            --tmp-dir ~/tmp/{rule} \
            --create-output-bam-index false"""

## Variant calling
# Modified input to skip recalibration
rule call_variants:
    """Call variants to make VCF file."""
    input:
        ref = config["ref_fasta"],
        bam = config["results"] + "alignments/recalibrated/{sample}.bam",
        bam_idx = config["results"] + "alignments/recalibrated/{sample}.bam.bai",
    output:
        vcf = config["results"] + "gvcf/{sample}.g.vcf.gz",
        tbi = config["results"] + "gvcf/{sample}.g.vcf.gz.tbi",
    conda: "../envs/gatk.yaml"
    threads: 9  # 4 is default
    resources: nodes = 9
    shell: """
        gatk --java-options '-Xmx16g' HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            -ERC GVCF \
            --native-pair-hmm-threads {threads} \
            --tmp-dir ~/tmp/{rule}"""

# Consolidation of GVCFs
rule create_sample_map:
    """Create sample map that contains names and paths to all VCFs to be used in consolidate rule.
    Note: This output file will need to be deleted if changing what will be added in the consolidate rule."""
    input:
        gvcfs = expand("{results}gvcf/{sample}.g.vcf.gz", results=config["results"], sample=SAMPLE_NAMES),
    output:
        sample_map = config["results"] + "db/{workspace}.sample_map",
    params:
        in_path = config["results"] + "gvcf/",
        samples = SAMPLE_NAMES,
    threads: 1
    resources: nodes = 1
    run:
        with open(output.sample_map, "w") as sample_map:
            print("Type", type(input.gvcfs))
            for gvcf in input.gvcfs:
                print("gvcf", gvcf)
                sample = gvcf.split("/")[-1].split(".")[0]
                sample_map.write(f"{sample}\t{gvcf}\n")

rule find_chromosomes:
    """Find all chromosomes from reference genome. Each line of the output file is just one chromosome's name."""
    input:
        config["ref_fasta"],  # Or something different
    output:
        config["results"] + "db/chromosomes.list",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    # Still need to generalize for whatever the count of chromosomes are.
    shell: """
        echo {{1..20}} | awk -v RS=' ' '{{print $1}} END {{print "X\\nY\\nMT"}}' > {output}
        """

rule consolidate:
    """Combine the chromosomes of .g.vcf files into GenomicsDB datastore."""
    input:
        contigs = config["results"] + "db/chromosomes.list",
        sample_map = config["results"] + "db/{workspace}.sample_map",
    output:
        config["results"] + "db/created_{workspace}.txt",
    params:
        db = config["results"] + "db/{workspace}",
        parallel_intervals = 4,  # Higher value requires more memory and number of file descriptor able to be used at the same time
    threads: 2  # Just for opening multiple .vcf files at once.
    resources: nodes = 2
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
    input:
        ref = config["ref_fasta"],
        # Note: The actual text file isn't what is required, but the datastore directory.
        # This .txt file is, however, is created only after the datastore has finished being built.
        db = config["results"] + "db/created_{workspace}.txt",
    output:
        vcf = config["results"] + "joint_call/{workspace}.vcf.gz",
        vcf_idx = config["results"] + "joint_call/{workspace}.vcf.gz.tbi",
    params:
        db = config["results"] + "db/{workspace}",
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    # Versions of GATK starting with 4.2.3.0 stop using "./." for missing genotypes and instead use "0/0".
    # This difference effects downstream analysis when converting into PLINK format.
    # And so, an earlier version must be used.
    shell: """
        gatk-4.2.2.0 --java-options '-Xmx8g' GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{params.db} \
            -O {output.vcf}"""

rule joint_call_cohort_per_chromosome:
    """Use GenomicsDB to jointly call a VCF file."""
    input:
        ref = config["ref_fasta"],
        # Note: The actual text file isn't what is required, but the datastore directory.
        # This .txt file is, however, is created only after the datastore has finished being built.
        db = config["results"] + "db/created_{workspace}.txt",
    output:
        vcf = config["results"] + "joint_call/chr/{workspace}.chr{chr}.vcf.gz",
        vcf_idx = config["results"] + "joint_call/chr/{workspace}.chr{chr}.vcf.gz.tbi",
    params:
        db = config["results"] + "db/{workspace}",
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    shell: """
        gatk-4.2.2.0 --java-options '-Xmx16g' GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{params.db} \
            -O {output.vcf} \
            -L {wildcards.chr} \
        """

## Split SNPs and indels
rule subset_mode:
    """Split into SNP- or indel-only .vcf. The wildcard `mode` can be "SNP" or "indel". Then keeps only biallelic sites."""
    input:
        vcf = config["results"] + "joint_call/{workspace}.vcf.gz",
        vcf_index = config["results"] + "joint_call/{workspace}.vcf.gz.tbi",
        ref_fasta = config["ref_fasta"],
    output:
        split = config["results"] + "joint_call/chr/{workspace}.{mode}.biallelic.vcf.gz",
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
            -Oz \
            -o {output.split} \
        """

