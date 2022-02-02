"""Contain rules for counting the number of reads for each genes.
This file is designed for RNA-seq data."""

rule index_transcripts:
    """Index reference transcripts."""
    input: transcripts=config["resources"] + "rna_transcripts.fna.gz"
    output: index=directory(config["resources"] + "rna_transcripts.fna.gz.salmon")
    params: min_kmer=31  # Minimum acceptable length for a valid match
    threads: 2  # Default for `-p`
    # --`geneMap` is a file containing 
    shell: "salmon index \
                -t {input.transcripts} \
                -i {output.index} \
                -k {params.min_kmer} \
                -p {threads}"

rule count_genes:
    """Create a summary of gene counts for single-end reads. .map can be generated from a features_table file from NCBI."""
    input: index=config["resources"] + "rna_transcripts.fna.gz.salmon",
           read=lambda wildcards: sample_path("{sample}".format(sample=wildcards.sample)) + "{sample}.fastq.gz",
           map=config["resources"] + "transcript_to_gene.map"
           #gff=config["resources"] + "genomic.gff"
    output: quants=directory(config["results"] + "gene_counts/unfiltered_counts/{sample}.counts")
    threads: 8
    # "A" infers library type.
    shell: "salmon quant \
                -i {input.index} \
                --libType A \
                -r {input.read} \
                -o {output.quants} \
                -p {threads} \
                --geneMap {input.map}"

rule merge_counts:
    """Merge counts from different samples."""
    input: quants=expand("{results}gene_counts/{sample}.counts", results=config["results"], sample=SAMPLE_NAMES)
    output: merged_quants=config["results"] + "gene_counts/merged.counts"
    params: quants=lambda wildcards, input: "{" + ",".join(input.quants) + "}"
    threads: 1
    # `-gene` flag uses gene quantification instead of transcript
    shell: "salmon quantmerge \
                --quants {params.quants} \
                --column numreads \
                -o {output.merged_quants} \
                --genes"

rule only_gene_names:
    """Remove non-coding mRNAs and sort by gene."""
    input: quants=config["results"] + "gene_counts/unfiltered_counts/{sample}.counts"
    output: filtered_quants=config["results"] + "gene_counts/merged_filtered.counts"
    threads: 1
    shell: "awk 'NR == 1; NR > 1 && $1 !~ /^XR_/ && ($2!=0 && $2!=0 && $3!=0 && $4!=0) {print $0 | \"sort\"}' merged_filtered.counts > {output.filtered_quants}"
