"""Rules related to determining genomic coverage of sequences."""

CONFIG = config["coverage"]
DIR = config["results"] + "coverage/"

rule find_coverage:
    """Find coverage of sequence. Helpful for finding what regions are being sequenced."""
    input: bam = config["results"] + "alignments/{sample}.bam",
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

rule sum_coverages:
    """Sum counts for each sample. Then filter out reads below a cutoff count."""
    input: DIR + "merged_coverage.bg",
    output: DIR + "summed_coverage.bg",
    params: cutoff = CONFIG["cutoff"],
    threads: 1
    conda: "../envs/coverage.yaml"
    # Find union of of `.bg` files.
    # `awk` sums from the fourth column to the end of line.
    shell: "awk '{{for (i=4; i<NF; i++) j+=$i; print $1,$2,$3,j; j=0}}' {input} \
            | awk '$4 > {params.cutoff} {{print $0}}' \
            > {output}"
