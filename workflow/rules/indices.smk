"""Rules for creating various types of index files."""

rule index_ref:
    """Create BWA index files for reference genome."""
    input: config["ref_fasta"],
    output: multiext(config["ref_fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: "bwa index {input}"

rule index_ref_fai:
    """Create .fai index for reference genome. Requires bgzipped reference."""
    input: config["ref_fasta"],
    output: config["ref_fasta"] + ".fai",
    conda: "../envs/bio.yaml"
    shell: "samtools fqidx {input}"

rule create_ref_dict:
    """Create .dict file for compressed reference genome."""
    #wildcard_constraints: fasta = "fasta|fna|fa"
    input: fasta = config["ref_fasta"],
    output: dict = ".".join(config["ref_fasta"].split(".")[:-2]) + ".dict",  # Replaces the ".fna.gz" ending with ".dict"
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    shell: "gatk CreateSequenceDictionary \
            -R {input.fasta}"

rule tbi_index:
    """Create .tbi index."""
    input: "{path}/{name}.vcf.gz",
    output: "{path}/{name}.vcf.gz.tbi",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: "bcftools index {input} \
                --tbi"

rule csi_index:
    """Create .csi index."""
    input: "{path}/{name}.vcf.gz",
    output: "{path}/{name}.vcf.gz.csi",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: "bcftools index {input} \
                --csi"

rule bai_index:
    """Create .bai index for .bam file."""
    input: "{path}/{sample}.bam",
    output: "{path}/{sample}.bam.bai",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: "samtools index {input}"