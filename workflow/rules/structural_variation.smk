"""Workflow for determing structural variation in genomes, that is, deletions, insertions, inversions, tandem duplications, andBND."""


## Short-read SV calling

rule call_SVs:
    """Call SVs per sample."""
    input:
        bam = config["results"] + "alignments/recalibrated/{sample}.bam",
        ref_fasta = config["ref_fasta"],
    output:
        bcf = config["results"] + "structural_variants/SVs/{sample}.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly.yaml"
    shell: """
        delly call {input.bam} \
            -g {input.ref_fasta} \
        > {output.bcf} \
        """

rule merge_SVs:
    """Merge BCFs with structural variants."""
    input:
        bcfs = expand(config["results"] + "structural_variants/SVs/{sample}.bcf", sample=SAMPLE_NAMES),
    output:
        bcf = config["results"] + "structural_variants/SVs/merged/{dataset}.ungenotyped.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly.yaml"
    shell: """
        delly merge {input.bcfs} \
            -o {output.bcf} \
        """

rule genotype_merged_SVs:
    """Genotype merged BCFs with structural variants."""
    input:
        bam = config["results"] + "alignments/recalibrated/{sample}.bam",
        bcf = config["results"] + "structural_variants/SVs/merged/{dataset}.ungenotyped.bcf",
        ref_fasta = config["ref_fasta"],
    output:
        bcf = config["results"] + "structural_variants/SVs/all_sites/{dataset}.{sample}.genotyped.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly.yaml"
    shell: """
        delly call {input.bam} \
            -g {input.ref_fasta} \
            -v {input.bcf} \
            -o {output.bcf} \
        """

rule merge_all_sites_SVs:
    """Merge samples with sites from all other samples."""
    input:
        bcfs = lambda wildcards: expand(
            config["results"] + "structural_variants/SVs/all_sites/{dataset}.{sample}.genotyped.bcf",
                dataset=wildcards.dataset,
                sample=SAMPLE_NAMES),
    output:
        bcf = config["results"] + "structural_variants/SVs/merged/{dataset}.genotyped.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools merge {input.bcfs}\
            -m id \
            -Ob \
            -o {output.bcf} \
        """

rule filter_SVs:
    """Filter structural variants."""
    input:
        bcf = config["results"] + "structural_variants/SVs/merged/{dataset}.genotyped.bcf",
    output:
        bcf = config["results"] + "structural_variants/SVs/merged/{dataset}.genotyped.filtered.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly.yaml"
    shell: """
        delly filter {input.bcf} \
            -f germline \
            -o {output.bcf} \
        """

rule split_by_SV_type:
    """Split merged BCF by type of SV."""
    wildcard_constraints:
        SV_type = "DEL|INS|DUP|INV|BND|ALL",
    input:
        bcf = config["results"] + "structural_variants/SVs/merged/{dataset}.genotyped.filtered.bcf",
    output:
        bcf = config["results"] + "structural_variants/SVs/merged/{dataset}.{SV_type}.genotyped.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools view {input.bcf} \
            -i "ALT='<{wildcards.SV_type}>'" \
            -Ob \
            -o {output.bcf} \
        """

## Mappability map

rule chop_ref_fna:
    """Chop reference FASTA into R1 and R2 FASTQ."""
    input:
        ref_fasta = config["ref_fasta"],
    output:
        R1 = config["results"] + "structural_variants/mappability/chopped_ref_fna.R1.fq.gz",
        R2 = config["results"] + "structural_variants/mappability/chopped_ref_fna.R2.fq.gz",
    params:
        R1 = config["results"] + "structural_variants/mappability/chopped_ref_fna.R1.fq.gz",
        R2 = config["results"] + "structural_variants/mappability/chopped_ref_fna.R2.fq.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly.yaml"
    shell: """
        dicey chop {input.ref_fasta} \
            --fq1 {params.R1} \
            --fq2 {params.R2} \
        """

rule align_chopped_ref_fna:
    """Align chopped reference FASTA."""
    input:
        ref_fasta = config["ref_fasta"],
        ref_indices = multiext(config["ref_fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
        R1 = config["results"] + "structural_variants/mappability/chopped_ref_fna.R1.fq.gz",
        R2 = config["results"] + "structural_variants/mappability/chopped_ref_fna.R2.fq.gz",
    output:
        bam = config["results"] + "structural_variants/mappability/self_aligned.bam",
    params:
        bwa_threads = 30,
        samtools_threads = 4,
    threads: 34
    resources: nodes = 34
    conda: "../envs/bio.yaml"
    shell: """
        bwa mem {input.ref_fasta} {input.R1} {input.R2} \
            -t {params.bwa_threads} \
        | samtools sort \
            -@ {params.samtools_threads} \
            -o {output.bam} \
        """

rule create_mappability_map:
    """Create mappability map for calling CNVs."""
    input:
        bam = config["results"] + "structural_variants/mappability/self_aligned.bam",
        bai = config["results"] + "structural_variants/mappability/self_aligned.bam.bai",
    output:
        map = config["results"] + "structural_variants/mappability/map.fa.gz",
    params:
        map = config["results"] + "structural_variants/mappability/map.fa",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly.yaml"
    shell: """
        dicey mappability2 {input.bam} \
            -o {output.map}; \
        gunzip {output.map} && bgzip {params.map} \
        """

## Copy number variation calling

rule call_CNVs:
    """Call CNVs per sample."""
    input:
        bam = config["results"] + "alignments/recalibrated/{sample}.bam",
        ref_fasta = config["ref_fasta"],
        map = config["results"] + "structural_variants/mappability/map.fa.gz",
        fai = config["results"] + "structural_variants/mappability/map.fa.gz.fai",
    output:
        bcf = config["results"] + "structural_variants/CNVs/{sample}.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly.yaml"
    shell: """
        delly cnv {input.bam} \
            -g {input.ref_fasta} \
            -l {input.map} \
            -o {output.bcf} \
        """

rule merge_CNVs:
    """Merge CNVs."""
    input:
        bcfs = expand(config["results"] + "structural_variants/CNVs/{sample}.bcf", sample=SAMPLE_NAMES),
    output:
        bcf = config["results"] + "structural_variants/CNVs/{dataset}.ungenotyped.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly.yaml"
    shell: """
        delly merge {input.bcfs} \
            --cnvmode \
            --pass \
            --minsize 1000 \
            --maxsize 100000 \
            -o {output.bcf} \
        """

rule genotype_merged_CNVs:
    """Genotype merged BCFs with copy number variants."""
    input:
        bam = config["results"] + "alignments/recalibrated/{sample}.bam",
        bcf = config["results"] + "structural_variants/merged/{dataset}.ungenotyped.bcf",
        ref_fasta = config["ref_fasta"],
        map = config["results"] + "structural_variants/mappability/map.fa.gz",
        fai = config["results"] + "structural_variants/mappability/map.fa.gz.fai",
    output:
        bcf = config["results"] + "structural_variants/CNVs/all_sites/{dataset}.{sample}.genotyped.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly.yaml"
    shell: """
        delly cnv {input.bam} \
            --segmentation \
            -v {input.bcf} \
            -g {input.ref_fasta} \
            -m {input.map} \
            -o {output.bcf} \
        """

rule merge_all_sites_CNVs:
    """Merge samples with sites from all other samples."""
    input:
        bcfs = lambda wildcards: expand(
            config["results"] + "structural_variants/CNVs/all_sites/{dataset}.{sample}.genotyped.bcf",
                dataset=wildcards.dataset,
                sample=SAMPLE_NAMES),
    output:
        bcf = config["results"] + "structural_variants/CNVs/merged/{dataset}.genotyped.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools merge {input.bcfs} \
            -m id \
            -Ob \
            -o {output.bcf} \
        """

rule filter_CNVs:
    input:
        bcf = config["results"] + "structural_variants/CNVs/merged/{dataset}.genotyped.bcf",
    output:
        bcf = config["results"] + "structural_variants/CNs/merged/{dataset}.genotyped.filtered.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly.yaml"
    shell: """
        delly classify {input.bcf} \
            --filter germline \
            -o {output.bcf} \
        """
