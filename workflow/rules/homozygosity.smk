rule count_length_of_contigs:
    """Works with NCBI FASTQ headers."""
    input:
        ref_fna = config["ref_fasta"],
    output:
        contig_lengths = config["results"] + "relatedness/roh/contig_lengths.tsv",
    shell: """
        gunzip -c {input.ref_fna} \
        | grep ">" \
        | cut -d ":" -f 4,6 \
        | sed 's/:/\t/' \
        > {output.contig_lengths} \
        """

rule scikit_allel_ROH:
    """Calculate runs of homozygosity. Reading pickle file from different version of Python than what it was created in can give a protocol error."""
    input:
        vcf = config["results"] + "haplotypes/SHAPEIT4/with_scaffold/{dataset}.{mode}.chr{chr}.vcf.gz",
        contig_lengths = config["results"] + "relatedness/roh/contig_lengths.tsv",
    output:
        roh_pickle = config["results"] + "relatedness/roh/{dataset}.{mode}.chr{chr}.roh_poisson.pickle",  # Stores a pandas datafrane
        froh_pickle = config["results"] + "relatedness/roh/{dataset}.{mode}.chr{chr}.froh_poisson.pickle",  # Stores a dictionary
    threads: 1
    resources: nodes = 1
    conda: "../envs/scikit.yaml"
    #script: "../scripts/runs_of_homozygosity.py"
    shell: """
        python3 workflow/scripts/runs_of_homozygosity.py {input.vcf} {input.contig_lengths} {output.roh_pickle} {output.froh_pickle} {wildcards.chr} \
        """