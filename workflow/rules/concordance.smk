"""Rules related to determining concordance between samples."""


# rule separate_sample:
#     """Create single-sample VCF file."""
#     input:
#         #joint_vcf = config["results"] + "genotypes/filtered/{dataset}_labels.{mode}.chr{chr}.filtered.vcf.gz",
#         joint_vcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
#     output:
#         indiv_vcf = config["results"] + "concordance/indiv_vcfs/{dataset}.{sample}.{mode}.chr{chr}.vcf.gz",
#     shell: """
#         bcftools view {input.joint_vcf} \
#             -s {wildcards.sample} \
#             -Oz \
#             -o {output.indiv_vcf} \
#         """

# rule merge_two_samples_at_common_sites:
#     """Merge samples from different VCF files."""
#     input:
#         vcf1 = config["results_main"] + "concordance/indiv_vcfs/GBS_WES_WGS_labels.{sample1}.{mode}.chr{chr}.vcf.gz",
#         tbi1 = config["results_main"] + "concordance/indiv_vcfs/GBS_WES_WGS_labels.{sample1}.{mode}.chr{chr}.vcf.gz.tbi",
#         vcf2 = config["results"] + "concordance/indiv_vcfs/{dataset}.{sample2}.{mode}.chr{chr}.vcf.gz",
#         tbi2 = config["results"] + "concordance/indiv_vcfs/{dataset}.{sample2}.{mode}.chr{chr}.vcf.gz.tbi",
#     output:
#         vcf = config["results"] + "concordance/merged_vcfs/{dataset}.{sample1}_{sample2}.{mode}.chr{chr}.vcf.gz",
#     shell: """
#         bcftools merge {input.vcf1} {input.vcf2} \
#             -Ou \
#         | bcftools view \
#             -i "F_MISSING=0" \
#             --min-ac 1 \
#             -Oz \
#             -o {output.vcf} \
#         """

# rule genotype_concordance:
#     """Calculate genotype concordance between two samples."""
#     input:
#         #vcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
#         vcf = config["results"] + "concordance/merged_vcfs/{dataset}.{sample1}_{sample2}.{mode}.chr{chr}.vcf.gz",
#     output:
#         concordance = config["results"] + "concordance/matches/{dataset}.{sample1}_{sample2}.{mode}.chr{chr}.concordance.tsv",
#     threads: 1
#     resources: nodes = 1
#     conda: "../envs/common.yaml"
#     # Reduce columns
#     # Remove header
#     # Keep only GT columns
#     # Replace "|" with "/"
#     # Swap any "1/0" with "0/1"
#     # Count matches and mismatches
#     shell: """
#         bcftools annotate {input.vcf} \
#             -x INFO,FORMAT \
#         | bcftools view \
#             -s {wildcards.sample1},{wildcards.sample2} \
#         | bcftools view \
#             -i 'F_MISSING=0' -H \
#         | cut -f 10-11 \
#         | sed 's/|/\//g' \
#         | sed 's/1\/0/0\/1/g' \
#         | awk -v OFS='\t' \
#             'BEGIN {{print "SAMPLE1","SAMPLE2","MATCHES","MISMATCHES","CONCORDANCE"}} \
#             {{if ($1 == $2) MATCH+=1; else MISMATCH+=1}} \
#             END {{CONCORDANCE = MATCH/(MATCH + MISMATCH); print "{wildcards.sample1}","{wildcards.sample2}",MATCH,MISMATCH,CONCORDANCE}}' \
#         > {output.concordance} \
#         """

# rule combine_concordances:
#     """Combine concordances information from different samples."""
#     input:
#         lambda wildcards: expand(config["results"] + "concordance/matches/{dataset}.WGS{animal_id}_GBS{animal_id}.{mode}.chr{chr}.concordance.tsv",
#             dataset=wildcards.dataset,
#             animal_id=SAMPLES,
#             mode=wildcards.mode,
#             chr=wildcards.chr),
#     output: config["results"] + "concordance/{dataset}.{mode}.chr{chr}.concordance.tsv",
#     # Take header
#     # Take each file, remove header, and add 0's to empty field
#     shell: """
#             cat {input[0]} | head -n 1 > {output}; \
#             for FILE in {input}; do \
#                 cat $FILE | sed '1d' | sed 's/\t\t/\t0\t/g' >> {output}; \
#             done \
#         """

# rule count_non_missing:
#     """Generate a table of counts."""
#     input:
#         vcf = config["results"] + "/genotypes/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
#     output:
#         called_count = config["results"] + "concordance/counts/{sample}.txt",
#     shell: """
#         bcftools view {input.vcf} \
#             -s {wildcards.sample} \
#         | bcftools view \
#             -i "F_MISSING=0" -H \
#         | wc -l \
#         > {output.called_count} \
#         """

# rule tabulate_non_missing:
#     """Generate a table of counts."""
#     input:
#         called_counts = expand(config["results"] + "concordance/counts/{sample}.txt", sample=SAMPLES),
#     output:
#         combined_counts = config["results"] + "concordance/counts.list",
#         tsv = config["results"] + "concordance/counts.tsv",
#     shell: """
#         cat {input.called_counts} > {output.combined_counts}; \
#         paste {config[samples]} {output.combined_counts} > {output.tsv} \
#         """

# ## Creating delimited file showing annotations for matching an nonmatching
# ## to be used to generate graphs

# rule extract_data:
#     input:
#         vcf = config["results"] + "GBS_with_markdup/concordance/merged_vcfs/GBS_with_markdup.WGS{organism_id}_GBS{organism_id}.SNP.chr{chr}.vcf.gz",
#     output:
#     run:
#         callset = allel.read_vcf(input.vcf,
#             ['variants/CHROM', 'variants/POS', 'variants/QUAL', 'variants/QD', 'variants/SOR', 'variants/FS', 'variants/MQ', 'variants/MQRankSum', 'variants/ReadPosRankSum', 'samples', 'calldata/GT'])



## For `bcftools gtcheck`


# Need to fix this now that `config["samples"]` doesn't exist anymore
# rule seq_pairs:
#     """Create TSV file match WGS to WES of the same id."""
#     input:
#         samples = config["samples"],
#     output:
#         pairs = config["results"] + "concordance/gtcheck/{dataset}.pairs.list",
#     run:
#         import polars as pl

#         pl.read_csv(input.samples, has_header=False, new_columns=["sample"]).with_columns(
#             pl.col("sample").str.slice(0, 3).alias("seq_type"),
#             pl.col("sample").str.split("S").list.last().str.split(".").list.first().alias("id"),
#         ).groupby("id").agg(pl.exclude("id")
#         ).select("sample").filter(
#             pl.col("sample").list.lengths() == 2
#         ).select(
#             pl.col("sample").list.get(0).alias("sample1"),
#             pl.col("sample").list.get(1).alias("sample2"),
#         ).sort(
#             "sample1", "sample2"
#         ).write_csv(output.pairs, has_header=False, separator="\t")


#######
# RTG #
#######

rule rtg_template:
    """Create SDF file, which is the required format for `rtg`."""
    input:
        ref_fasta = config["ref_fasta"],
    output:
        sdf = config["results"] + "concordance/rtg/sdf/summary.txt",
    params:
        sdf = config["results"] + "concordance/rtg/sdf",
    conda: "../envs/rtg.yaml"
    shell: """
        rtg format {input.ref_fasta} \
            --format fasta \
            -o {params.sdf} \
        """

# This rule gives an error about the output directory already existing. Might be from Snakemake generating an output file when it starts the run
rule rtg_vcfeval:
    """Compare a pair of samples between VCFs. Use LRS as "baseline" and WGS as "calls" to compare against it."""
    input:
        baseline = config["results"] + "structural_variants/delly/merged/LRS.vcf",
        calls = config["results"] + "structural_variants/delly/merged/U42_WGS.vcf",
        sdf = config["results"] + "concordance/rtg/sdf/summary.txt",
    output:
        config["results"] + "concordance/rtg/vcfeval/{sample}/{sample}.txt",  # To change
    params:
        out_dir = config["results"] + "concordance/rtg/vcfeval/{sample}",
        sdf_dir = config["results"] + "concordance/rtg/sdf",
    conda: "../envs/rtg.yaml"
    shell: """
        rtg vcfeval \
            --baseline {input.baseline} \
            --calls {input.calls} \
            -o {params.out_dir} \
            --sample LRS{wildcards.sample},WGS{wildcards.sample} \
            --template {params.sdf_dir} \
        """

##########
# SV-Pop #
##########

# rule sv_pop:
#     """"""




# rule bcftools_stats:
#     """Find concordance."""

rule gtcheck:
    """Find concordance between samples of same individual in the same VCF."""
    wildcard_constraints:
        method = "abstract|discrete",
        coverage = "genome|exome",
    input:
        # vcfs = lambda wildcards: expand(config["results"] + "hard_filtered/pass/{dataset}.SNP.chr{chr}.vcf.gz",
        #     dataset=wildcards.dataset,
        #     chr=[str(i) for i in range(1,21)]
        # ),
        # tbis = lambda wildcards: expand(config["results"] + "hard_filtered/pass/{dataset}.SNP.chr{chr}.vcf.gz.tbi",
        #     dataset=wildcards.dataset,
        #     chr=[str(i) for i in range(1,21)]
        # ),
        vcf = lambda wildcards: expand(config["results"] + "hard_filtered/pass/{dataset}.SNP.autosomal.vcf.gz",
            dataset=wildcards.dataset,

        ),
        tbi = lambda wildcards: expand(config["results"] + "hard_filtered/pass/{dataset}.SNP.autosomal.vcf.gz.tbi",
            dataset=wildcards.dataset,
        ),
        exons = config["resources"] + "annotations/exons.tsv"
    output:
        discord = config["results"] + "concordance/gtcheck_{method}_{coverage}/{dataset}.{individual}.discord.txt",
    params:
        samples = lambda wildcards: collect_sample_libraries(wildcards.individual),
        add_options = lambda wildcards: "" if wildcards.method == "abstract" else "-u GT,GT -E 0",
        optional_exome = lambda wildcards: "" if wildcards.coverage == "genome" else "-T {input.exons}",
    #-P {input.pairs} \
    #            -u GT,GT \
    #        -e 0 \
            # bcftools concat {input.vcfs} \
            # -n \
            # -Ou \
    shell: """
        bcftools +setGT {input.vcf} \
            -Ou \
            -- \
            -tq \
            -i 'FORMAT/GQ<30' \
            -n . \
        | bcftools view \
            -e 'FORMAT/DP=0' \
            {params.optional_exome} \
            -s {params.samples} \
            -Ou \
        | bcftools gtcheck \
            {params.add_options} \
        > {output.discord} \
        """



# rule find_replicate_ids:
#     input:
#         runs = config["runs"]
#     output:
#         replicates = temp(config["resources"] + "samples/replicates.list")
#     # Only works if the batch name does not have an underscore
#     shell: """
#         cat {input.runs} \
#         | cut -d '_' -f 1,2 \
#         | sort \
#         | uniq \
#         | cut -d '/' -f 2 \
#         | cut -d _ -f 1 \
#         | cut -c 4- \
#         | sort \
#         | uniq -d \
#         > {output.replicates} \
#         """

# with open(config["runs"]) as f:
#     SAMPLE_RUNS = f.read().splitlines()
import collections
animals_with_dup = [sample.split("_")[0][3:] for sample in SAMPLES]
REPLICATE_IDS = [item for item, count in collections.Counter(animals_with_dup).items() if count > 1]

rule merge_gtcheck:
    """Merge results from gtcheck."""
    wildcard_constraints:
        method = "abstract|discrete",
        coverage = "genome|exome",
    input:
        gtchecks = lambda wildcards: expand(config["results"] + "concordance/gtcheck_{method}_{coverage}/{dataset}.{individual}.discord.txt",
            method=wildcards.method,
            coverage=wildcards.coverage,
            dataset=wildcards.dataset,
            #individual=['30009','16653','16356','18401','17535','17553','32351','16652','17730'],
            individual=REPLICATE_IDS,
        )
    output:
        merged = config["results"] + "concordance/gtcheck_{method}_{coverage}/{dataset}.discord.txt",
    shell: """
        cat {input.gtchecks} \
        | grep ^DC \
        > {output.merged} \
        """
# Fix to include header

# rule gtcheck_imputed_WGS_vs_WES:
#     """Find concordance."""
#     input:
#         vcf = config["results"] + "haplotypes/SHAPEIT5_WES/{dataset}.{mode}.autosomal.bcf",
#         tbi = config["results"] + "haplotypes/SHAPEIT5_WES/{dataset}.{mode}.autosomal.bcf.csi",
#         genotypes = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.{mode}.autosomal.bcf",
#         genot_tbi = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.{mode}.autosomal.bcf.csi",
#         pairs = config["results"] + "concordance/gtcheck/{dataset}.pairs.list",
#     output:
#         discord = config["results"] + "concordance/gtcheck/{tool}/{dataset}.SNP.discord.txt",
#     shell: """
#         bcftools gtcheck {input.vcf} \
#             -u GT,GT \
#             -g {input.genotypes} \
#             -e 0 \
#             -P {input.pairs} \
#         > {output.discord} \
#         """

rule gtcheck_percent:
    """Find percent discordance from gtcheck output."""
    input:
        discord = config["results"] + "concordance/gtcheck/{tool}/{dataset}.{mode}.discord.txt",
    output:
        discord = config["results"] + "concordance/gtcheck/{tool}/{dataset}.{mode}.discord_percent.txt",
    shell: """
        grep ^DC {input.discord} \
        | awk '{{print $1,$2,$3,$4/$6}}' \
        > {output.discord} \
        """

## For vcftools-based concordance

rule vcftools_genotype_counts:
    """Create matrix detailing count of reference haplotypes, that is 0/0 -> 0, 0/1 -> 1, 1/1 -> 2."""
    input:
        vcf = config["results"] + "haplotypes/SHAPEIT4/with_scaffold/{dataset}.{mode}.vcf.gz",
    output:
        matrix = config["results"] + "concordance/vcftools/SHAPEIT4/{dataset}.{mode}.012",
        indv = config["results"] + "concordance/vcftools/SHAPEIT4/{dataset}.{mode}.012.indv",
        pos = config["results"] + "concordance/vcftools/SHAPEIT4/{dataset}.{mode}.012.pos",
    params:
        out_prefix = lambda wildcards, output: ".".join(output.matrix.split(".")[:-1]),
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        vcftools \
            --gzvcf {input.vcf} \
            --012 \
            --out {params.out_prefix} \
        """

rule matrix_discordance:
    """Calculate discordance from vcftools matrix."""
    input:
        matrix = config["results"] + "concordance/vcftools/SHAPEIT4/{dataset}.{mode}.discord.txt",
        indv = config["results"] + "concordance/vcftools/SHAPEIT4/{dataset}.{mode}.012.indv",
        pos = config["results"] + "concordance/vcftools/SHAPEIT4/{dataset}.{mode}.012.pos",
        pairs = config["results"] + "concordance/gtcheck/{dataset}.pairs.list",
    output:
        discord = config["results"] + "concordance/vcftools/SHAPEIT4/{dataset}.{mode}.discord.tsv",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        while read -r SAMPLE1 SAMPLE2; do \
            echo $SAMPLE1 $SAMPLE2;
            LINE1=$(grep -n $SAMPLE1 {input.indv} | cut -d ":" -f 1)
            LINE2=$(grep -n $SAMPLE2 {input.indv} | cut -d ":" -f 1)
            
        done < {input.pairs} \
        """
