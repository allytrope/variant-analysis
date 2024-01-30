"""Rules for creating various types of index files."""

rule index_ref:
    """Create BWA index files for reference genome."""
    input: config["ref_fasta"],
    output: multiext(config["ref_fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: "bwa index {input}"

rule fai_index:
    """Create .fai index for .fna.gz file. Must be bgzipped."""
    wildcard_constraints: fna = "fasta|fna|fa"
    input: "{path}/{name}.{fna}.gz"
    output: "{path}/{name}.{fna}.gz" + ".fai",
    conda: "../envs/bio.yaml"
    shell: "samtools faidx {input}"

rule gzi_index:
    """Create .gzi index for .fna.gz file. Must be bgzipped."""
    wildcard_constraints: fna = "fasta|fna|fa"
    input: "{path}/{name}.{fna}.gz"
    output: "{path}/{name}.{fna}.gz" + ".gzi",
    conda: "../envs/bio.yaml"
    shell: """
        bgzip {input} \
            -r \
        """

# Need to implement for Whatshap
# rule pyfaidx_index:
#     """Create .fai index for .fna.gz file.
#     Creates a different index than the other rule,
#     which is needed for WhatsHap."""
#     wildcard_constraints: fna = "fasta|fna|fa"
#     input: "{path}/{name}.{fna}.gz"
#     output: "{path}/faidx/{name}.{fna}.gz" + ".fai",
#     conda: "../envs/bio.yaml"
#     shell: "faidx {input}"

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

rule minimap_index:
    """Create .mmi index of reference fasta for Minimap2."""
    input:
        ref_fasta = config["ref_fasta"],
    output:
        mmi = config["ref_fasta"] + ".mmi",
    shell: """
        minimap2 {ref.fasta} \
            -d {output.mmi} \
        """

rule bai_index:
    """Create .bai index for .bam file."""
    input: "{path}/{name}.bam",
    output: "{path}/{name}.bam.bai",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: "samtools index {input}"

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
    wildcard_constraints:
        ext = "bcf|vcf.gz",
    input: "{path}/{name}.{ext}",
    output: "{path}/{name}.{ext}.csi",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: "bcftools index {input} \
                --csi"
