"""Rules related to determining genomic coverage."""

CONFIG = config["coverage"]
DIR = config["results"] + "coverage/"

## Determining common loci
rule find_coverage:
    """Find coverage of sequence. Helpful for finding what regions are being sequenced."""
    input: bam = config["results"] + "alignments_recalibrated/{sample}.bam",
    output: counts = DIR + "per_read/{sample}.bg",
    threads: 1
    conda: "../envs/bio.yaml"
    shell: "bedtools genomecov \
                -ibam {input.bam} \
                -bg > {output.counts}"

rule merge_coverages:
    """Combine coverage from multiple sequences. Creates new column of counts for each sample."""
    input: expand("{dir}per_read/{sample}.bg", dir=DIR, sample=SAMPLE_NAMES),
    output: DIR + "merged_coverage.bg",
    threads: 1
    conda: "../envs/bio.yaml"
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
    conda: "../envs/bio.yaml"
    shell: """
        awk -v 'OFS=\t' '{{for (i=4; i<=NF; i++) {{if ($i > 0) count += 1}}; if (count/(NF - 3) > {params.cutoff}) print $1,$2,$3; count = 0}}' {input}
        | grep '^NC'
        | bedtools merge
            -d 1 
        > {output}"""

## Finding fraction intersecting reference .bed
rule intersect_bed_count:
    """Find the number of intersections of .bam and .bed file."""
    input:
        target = "/data/infectious/malaria/marmoset/align/{sample}/{sample}_resorted_realigned_rg_recal.bam",
        ref_bed = config["resources"] + "ref_fasta/C_jacchus3.2.1.cdna.all.bed",
    output: intersection = DIR + "intersections/{sample}.stats",
    threads: 1
    conda: "../envs/bio.yaml"
    shell: """
        bedtools intersect \
            -a {input.target} \
            -b {input.ref_bed} \
        | samtools view \
        | wc -l \
        > {output}; \
        samtools view {input.target} \
        | wc -l \
        >> {output} \
        """

rule intersectional_fraction:
    """Merge counts."""
    input:
        #total = "/data/infectious/malaria/marmoset/align/{sample}/{sample}_resorted_realigned_rg_recal.bam",
        intersections = expand(DIR + "intersections/{sample}.stats", DIR=DIR, sample=SAMPLE_NAMES),
    output: total = DIR + "intersections/total.stats",
    threads: 1
    conda: "../envs/bio.yaml"
    shell: """
        for i in {{1..2}}; \
        do awk -v i=$i 'FNR==i {{ print }}' {input.intersections} | awk '{{ sum+=$1 }} END {{ print sum }}'; \
        done \
        > {output.total} \
        """

## VCF-BED intersection
rule vcf_bed_intersection:
    """Keep variants in VCF only if in .bed file."""
    input: vcf = config["results"] + "genotypes/filtered/GBS_WES_WGS.SNP.chr{chr}.filtered.vcf.gz",
        bed = config["results"] + "coverage/common_loci.100.bed",
        subset = config["resources"] + "samples/WES_WGS_no_label.list",  # Just temporary
    output: vcf = config["results"] + "genotypes/filtered/WES_WGS.SNP.chr{chr}.WES_regions.min-ac_3.vcf.gz",
    threads: 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools view {input.vcf} \
            -R {input.bed} \
            -S {input.subset} \
            --min-ac 3 \
            -Oz \
            -o {output.vcf} \
        """

rule nonmissing_genotypes:
    """Find fraction of genotypes not missing in samples of .vcf file."""
    input: vcf = config["results"] + "coverage/{dataset}.test.vcf",  # Set this to whatever .vcf file is of interest.
    output: config["results"] + "coverage/{dataset}.GT_stats.csv",
    #conda: "../envs/bio.yaml",
    run:
        import allel
        import pandas as pd

        callset = allel.read_vcf(input.vcf, ['samples', 'calldata/GT'])
        g = allel.GenotypeArray(callset['calldata/GT'])
        nonmissing_GT = g.count_called(axis=0)/len(g)

        data = {
            "sample": callset['samples'],
            "nonmissing_GT": nonmissing_GT,
        }
        with open(output, 'w') as f:
            f.write("# Based on the file: " + output)
        df = pd.DataFrame(data).set_index('sample')
        df.to_csv(output, mode="a")

# rule proportions_of_matching_WGS_and_GBS:
#     """Create a table for finding the percent of sites that match between WGS and GBS for the same sample."""
#     input: vcf = config["results"] + "genotypes/filtered/GBS_WES_WGS_labels.SNP.chr11.filtered.vcf.gz",
#     output: vcf = "",
#     threads: 1
#     resources: nodes = 1
#     run:
#         from collections import Counter

#         import allel
#         import pandas as pd

#         callset = allel.read_vcf(input.vcf, ['samples', 'calldata/GT'])
#         g = allel.GenotypeArray(callset['calldata/GT'])

#         samples = list(callset['samples'])
#         animal_ids = None
#         for animal_id in animal_ids:
#             try:
#                 gbs_idx = samples.index(f'GBS{animal_id}')
#                 wgs_idx = samples.index(f'WGS{animal_id}')
#             except:
#                 continue
            
#             counter = Counter([(i, j) for i, j in g[:, gbs_idx] == g[:, wgs_idx]])
#             matching_rate = counter[(True, True)]/counter.total()




#             if f"GBS{animal_id} in "callset['samples'] and f"WGS{animal_id}" in callset['samples']:
#                 gbs_gt = g[:, ]

rule windowed_depth:
    """Find depth across windows in one sample."""
    input:
        bam = config["results"] + "alignments/recalibrated/{sample}.bam",
    output:
        global_dist = config["results"] + "coverage/mosdepth/{sample}.mosdepth.global.dist.txt",
        region_dist = config["results"] + "coverage/mosdepth/{sample}.mosdepth.region.dist.txt",
        summary = config["results"] + "coverage/mosdepth/{sample}.mosdepth.summary.txt",
        bed = config["results"] + "coverage/mosdepth/{sample}.regions.bed.gz",
        csi = config["results"] + "coverage/mosdepth/{sample}.regions.bed.gz.csi",
    params:
        prefix = lambda wildcards, output: output.bed.split(".")[0],
        window_size = 5_000_000,
    threads: 4  # Docs say to use 4 or fewer
    resources: nodes = 4
    conda: "../envs/bio.yaml"
    shell: """
        mosdepth {params.prefix} {input.bam} \
            --by {params.window_size} \
            -n \
            --fast-mode \
            -t {threads} \
        """

rule merge_windowed_depth:
    """Merge windowed depth into one file for viewing joint figure."""
    input:
        #bed = expand(config["results"] + "coverage/mosdepth/{sample}.regions.bed.gz", sample=SAMPLE_NAMES),
        bed = expand(config["results"] + "coverage/mosdepth/{sample}.regions.bed.gz",
            sample=SAMPLE_NAMES),
    output:
        merged_bed = config["results"] + "coverage/mosdepth/merged.bed",
    shell: """
        for BED in {input.bed}; do \
            SAMPLE=$(echo $BED | cut -d "." -f 1 | rev | cut -d "/" -f 1 | rev); \
            gunzip -c $BED \
            | awk -v SAMPLE=$SAMPLE 'BEGIN {{OFS="\\t"}} {{print SAMPLE, $0}}' \
            >> {output.merged_bed}; \
        done \
        """

## for figure --> notebooks/coverage.ipynb
