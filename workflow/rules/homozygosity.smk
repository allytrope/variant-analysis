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
    wildcard_constraints:
        #chr="[0-9]+",
        seq = "WGS|WES|merged",
    input:
        # bcf = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.{mode}.chr{chr}.bcf",
        # csi = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.{mode}.chr{chr}.bcf.csi",
        #bcf = config["results"] + "haplotypes/SHAPEIT5_{seq}/{dataset}.{mode}.chr{chr}.bcf",
        #csi = config["results"] + "haplotypes/SHAPEIT5_{seq}/{dataset}.{mode}.chr{chr}.bcf.csi",
        bcf = config["results"] + "genotypes/pass/{dataset}.common_between_founding_cohorts2.SNP.chr{chr}.bcf",
        csi = config["results"] + "genotypes/pass/{dataset}.common_between_founding_cohorts2.SNP.chr{chr}.bcf.csi",
        #csi = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.vcf.gz",
        #bed = config["results"] + "coverage/common_WES_0.5_loci.bed",
    output:
        #roh = config["results"] + "roh/bcftools/{dataset}.{mode}.chr{chr}.roh",
        roh = config["results"] + "roh/bcftools/{dataset}.common_between_founding_cohorts2.chr{chr}.both.roh",
    #-R {input.bed} \
    #-G 30 \
    params:
        #rec_rate = 0.00433  # For rhesus
        #rec_rate = 0.003  # From rhesus WES discordance
        rec_rate = .00000043
    shell: """
        bcftools view {input.bcf} \
            --min-ac 1 \
            -Ou \
        | bcftools roh \
            --rec-rate {params.rec_rate} \
            -Or \
            -o {output.roh} \
        """

rule bcftools_ROH_split_rows:
    """Split rows starting with ST and RG."""
    wildcard_constraints:
        seq = "WGS|WES|merged",
        tag="ST|RG",
    input:
        roh = config["results"] + "roh/bcftools/{dataset}.{subset}.chr{chr}.both.roh",
    output:
        roh = config["results"] + "roh/bcftools/{dataset}.{subset}.chr{chr}.{tag}.roh",
    shell: """
        grep '^[{wildcards.tag}|#]' {input.roh} > {output.roh} \
        """

## PLINK ROH
# rule PLINK_ROH:
#     input:
#         # It is suggested that a pruned input is used.
#         bcf = config["results"] + "genotypes/pruned/plink/{dataset}.{subset}.{mode}.chr{chr}.bcf",
#     shell: """
        
#         """

## scikit-allel ##

rule scikit_allel_ROH:
    """Calculate runs of homozygosity. Reading pickle file from different version of Python than what it was created in can give a protocol error."""
    input:
        #vcf = config["results"] + "haplotypes/SHAPEIT5_merged/{dataset}.{mode}.chr{chr}.vcf.gz",
        #tbi = config["results"] + "haplotypes/SHAPEIT5_merged/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
        vcf = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz.tbi",
        contig_lengths = config["results"] + "relatedness/roh/contig_lengths.tsv",
    output:
        # roh_pickle = config["results"] + "roh/scikit-allel/{dataset}.{mode}.chr{chr}.roh_poisson.pickle",  # Stores a pandas datafrane
        # froh_pickle = config["results"] + "roh/scikit-allel/{dataset}.{mode}.chr{chr}.froh_poisson.pickle",  # Stores a dictionary
        roh_pickle = config["results"] + "roh/scikit-allel/{dataset}.{subset}.{mode}.chr{chr}.roh_poisson.tsv",
        froh_pickle = config["results"] + "roh/scikit-allel/{dataset}.{subset}.{mode}.chr{chr}.froh_poisson.tsv",
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
# Filter out low quality GQs
# Remove "END="
# Select START, END, and GT columns
# Set any sites missing END to START
# Split into new columns on "\" and "|"
# Create only columns for length and zygosity
# Count up heterozygous and total sites

rule heterozygosity_from_gvcf:
    """Count up number of heterozygous sites from gVCF."""
    input:
        gvcf = config["results"] + "gvcf/{sample}.g.vcf.gz",
    output:
        het = config["results"] + "heterozygosity/gvcf_counts/{sample}.het",
    shell: """
        bcftools view {input.gvcf} \
            -e 'GQ<30' \
        | bcftools annotate \
            -x ^INFO/END,FORMAT \
        | bcftools view \
            -H \
        | sed 's/END=//' \
        | cut \
            -f 2,8,10 \
        | awk \
            'BEGIN {{FS="\t";OFS="\t"}} \
            {{if ($2 == ".") $2=$1; print $0}}' \
        | sed 's/|/\t/;s/\//\t/' \
        | awk \
            'BEGIN {{FS="\t";OFS="\t"}} \
            {{LEN=$2-$1+1; if ($3 == $4) ZYG="hom"; else ZYG="het"; print LEN,ZYG}}' \
        | awk \
            'BEGIN {{FS="\t";OFS="\t"; print "SAMPLE","HET","TOTAL"}} \
            {{if ($2 == "hom") HOM+=$1; else if ($2 == "het") HET+=$1; TOTAL+=$1}} \
            END {{print "{wildcards.sample}",HET,TOTAL}}' \
        > {output}
        """

rule spawn_heterozygosity_jobs:
    """Spawn jobs for the rule `heterozygosity_from_gvcf`."""
    input:
        het = expand(config["results"] + "heterozygosity/gvcf_counts/{sample}.het",
            sample=collect_samples(fmt="{seq}{indiv}_{library}")
        )
    output:
        config["results"] + "heterozygosity/gvcf_counts.het",
    shell: """
        echo "SAMPLE\tHET\tTOTAL" > {output};
        cat {input.het} | grep -v '^SAMPLE' >> {output}
        """
