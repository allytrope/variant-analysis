"""Rules related to determining genomic coverage."""

CONFIG = config["coverage"]
DIR = config["results"] + "coverage/"

rule find_coverage:
    """Find coverage of sequence. Helpful for finding what regions are being sequenced."""
    input: bam = config["results"] + "alignments_recalibrated/{sample}.bam",
    output: counts = DIR + "per_read/{sample}.bg",
    threads: 1
    conda: "../envs/coverage.yaml"
    shell: "bedtools genomecov \
                -ibam {input.bam} \
                -bg > {output.counts}"

rule merge_coverages:
    """Combine coverage from multiple sequences. Creates new column of counts for each sample."""
    input: expand("{dir}per_read/{sample}.bg", dir=DIR, sample=SAMPLE_NAMES),
    output: DIR + "merged_coverage.bg",
    threads: 1
    conda: "../envs/coverage.yaml"
    shell: "bedtools unionbedg \
                -i {input} \
                > {output}"

# To be modified when using ENSEMBL reference (i.e., numbered chromosomes).
rule find_common_loci:
    """Find loci common to most samples based on cutoff value.
    This works by finding the number of samples with reads at that position, keeping only those above a given cutoff, and then merging intervals that overlap or are immediately adjacent."""
    input: DIR + "merged_coverage.bg",
    output: DIR + "common_loci.bed",
    params: cutoff = CONFIG["cutoff"],  # 0 <= cutoff <= 1
    threads: 1
    conda: "../envs/coverage.yaml"
    shell: """
        awk -v 'OFS=\t' '{{for (i=4; i<=NF; i++) {{if ($i > 0) count += 1}}; if (count/(NF - 3) > {params.cutoff}) print $1,$2,$3; count = 0}}' {input}
        | grep '^NC'
        | bedtools merge
            -d 1 
        > {output}"""
