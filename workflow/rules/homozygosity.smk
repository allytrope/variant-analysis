rule count_length_of_contigs:
    """Works with ENSEMBL FASTQ headers."""
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
        """\
        
## bcftools roh ###

rule bcftools_ROH:
    """Find ROHs using BCFtools."""
    input:
        bcf = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.{mode}.chr{chr}.bcf",
        csi = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.{mode}.chr{chr}.bcf.csi",
        bed = config["results"] + "coverage/common_WES_0.5_loci.bed",
    output:
        roh = config["results"] + "roh/SHAPEIT5_WGS/{dataset}.{mode}.chr{chr}.both.roh",
    #-R {input.bed} \
    shell: """
        bcftools roh {input.bcf} \
            -G 30 \
        > {output.roh} \
        """

rule bcftools_ROH_split_rows:
    """Split rows starting with ST and RG."""
    wildcard_constraints:
        tag="ST|RG",
    input:
        roh = config["results"] + "roh/SHAPEIT5_WGS/{dataset}.{mode}.chr{chr}.both.roh",
    output:
        roh = config["results"] + "roh/SHAPEIT5_WGS/{dataset}.{mode}.chr{chr}.{tag}.roh",
    shell: """
        grep '^[{wildcards.tag}|#]' {input.roh} > {output.roh} \
        """

## scikit-allel ##

rule scikit_allel_ROH:
    """Calculate runs of homozygosity. Reading pickle file from different version of Python than what it was created in can give a protocol error."""
    input:
        vcf = config["results"] + "haplotypes/SHAPEIT5_merged/{dataset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "haplotypes/SHAPEIT5_merged/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
        contig_lengths = config["results"] + "relatedness/roh/contig_lengths.tsv",
    output:
        roh_pickle = config["results"] + "relatedness/roh/merged/{dataset}.{mode}.chr{chr}.roh_poisson.pickle",  # Stores a pandas datafrane
        froh_pickle = config["results"] + "relatedness/roh/merged/{dataset}.{mode}.chr{chr}.froh_poisson.pickle",  # Stores a dictionary
    threads: 1
    resources: nodes = 1
    conda: "../envs/scikit.yaml"
    #script: "../scripts/runs_of_homozygosity.py"
    shell: """
        python3 workflow/scripts/runs_of_homozygosity.py {input.vcf} {input.contig_lengths} {output.roh_pickle} {output.froh_pickle} {wildcards.chr} \
        """
# For graphing .pickle files ---> notebooks/runs_of_homozygosity.ipynb


rule scikit_allel_diversity:
    """Estimate nucleotide diversity in windows. Output files are formatted to be graphed with the R tool chromoMap."""
    input: 
        vcf = config["results"] + "haplotypes/SHAPEIT4/with_scaffold/{dataset}.{mode}.chr{chr}.vcf.gz",
    output:
        annotations = config["results"] + "relatedness/diversity/diversity.tsv",
        chromosomes = config["results"] + "relatedness/diversity/chromosomes.tsv",
    conda: "../envs/scikit.yaml"
    #script: "../scripts/diversity.py"
    shell: """
        python3 workflow/scripts/diversity.py {input.vcf} {output.annotations} {output.chromosomes} \
    """
# For graphing files ---> notebooks/nucleotide_diversity.ipynb