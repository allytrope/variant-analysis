"""Rules related to determining genomic coverage."""

#CUTOFF = 0.8  # config["coverage"]
DIR = config["results"] + "coverage/"

## Determining common loci
rule find_coverage:
    """Find coverage of sequence. Helpful for finding what regions are being sequenced."""
    input:
        bam = config["results"] + "alignments/recalibrated/{batch}/{sample_run}.bam",
    output:
        counts = DIR + "per_read/{batch}/{sample_run}.bg",
    threads: 1
    conda: "../envs/bio.yaml"
    shell: """
        bedtools genomecov \
            -ibam {input.bam} \
            -bg > {output.counts}
        """

rule merge_coverages:
    """Combine coverage from multiple sequences. Creates new column of counts for each sample."""
    input:
        expand("{dir}per_read/{run_path}.bg", dir=DIR, run_path=SAMPLE_RUNS),
    output:
        DIR + "merged_coverage.bg",
        #DIR + "merged_WES_coverage.bg",
    threads: 1
    conda: "../envs/bio.yaml"
    shell: """
        bedtools unionbedg \
            -i {input} \
            > {output} \
        """

# To be modified when using ENSEMBL reference (i.e., numbered chromosomes).
rule find_common_loci:
    """Find loci common to most samples based on cutoff value.
    This works by finding the number of samples with reads at that position, keeping only those above a given cutoff, and then merging intervals that overlap or are immediately adjacent."""
    input:
        #DIR + "merged_coverage.bg",
        DIR + "merged_WES_coverage.bg",
    output:
        #DIR + "common_loci.bed",
        DIR + "common_WES_loci.bed",
    params: cutoff = 0.5,  # 0 <= cutoff <= 1
    threads: 1
    conda: "../envs/bio.yaml"
    shell: """
        awk -v 'OFS=\t' '{{ {{for (i=4; i<=NF; i++) {{if ($i > 0) count += 1}} }}; if (count/(NF - 3) > {params.cutoff}) print $1,$2,$3; count = 0}}' {input} \
        | bedtools merge \
            -d 1 \
        > {output} \
        """

# rule count_exon_length:
#     output: "total_exon_length.txt"
#     shell: """
#         cat GCF_009663435.1_Callithrix_jacchus_cj1700_1.1_genomic.chr_renamed.exons.no_header.gff | awk 'BEGIN {FS="\t"} {sum+=$5-$4+1} END {print sum}' |
        
#         """

rule intersect_exon_count:
    """Find count of bases where the .bam intersects with .bed file."""
    input:
        #target = "/data/infectious/malaria/marmoset/align/{sample}/{sample}_resorted_realigned_rg_recal.bam",
        target = config["results"] + "alignments/recalibrated/{batch}/{seq}{sample_run}.bam",
        #ref_bed = config["resources"] + "ref_fasta/C_jacchus3.2.1.cdna.all.bed",
        ref_bed = config["resources"] + "annotations/GCF_009663435.1_Callithrix_jacchus_cj1700_1.1_genomic.chr_renamed.exons.no_header.gff",
    output: intersection = DIR + "intersections/exons/{batch}/{seq}{sample_run}.stats",
    threads: 1
    conda: "../envs/bio.yaml"
    shell: """
        bedtools intersect \
            -a {input.target} \
            -b {input.ref_bed} \
        | samtools fasta \
        | grep -v '^>' \
        | tr -d '\n' \
        | wc -m \
        > {output.intersection} \
        """

rule combine_exon_counts:
    input:
        intersection = expand(DIR + "intersections/exons/{sample_path}.stats",
            DIR=DIR,
            sample_path=SAMPLE_RUNS,),
        total_length = DIR + "intersections/exons/total_exon_length.txt",
    output:
        file_names = DIR + "intersections/exons/file_names.txt",
        file_contents = DIR + "intersections/exons/file_contents.txt",
        total_exome_length = DIR + "intersections/exons/total_exome_length.txt",
    shell: """
        TOTAL_LENGTH=$(cat {input.total_length}); \
        for FILE in {input.intersection}; do \
            echo $FILE >> {output.file_names}; \
            cat $FILE >> {output.file_contents}; \
            cat $TOTAL_LENGTH >> {output.total_exome_length}; \
        done \
        """




## Finding fraction intersecting reference .bed
rule intersect_bed_count:
    """Find the number of intersections of .bam and .bed file."""
    input:
        #target = "/data/infectious/malaria/marmoset/align/{sample}/{sample}_resorted_realigned_rg_recal.bam",
        target = config["results"] + "alignments/recalibrated/{batch}/{sample_run}.bam",
        #ref_bed = config["resources"] + "ref_fasta/C_jacchus3.2.1.cdna.all.bed",
        ref_bed = config["resources"] + "annotations/GCF_009663435.1_Callithrix_jacchus_cj1700_1.1_genomic.chr_renamed.exons.gff"
    output: intersection = DIR + "intersections/{batch}/{sample}.stats",
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
        intersections = expand(DIR + "intersections/{sample}.stats", DIR=DIR, sample=SAMPLES),
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

## Depth
rule windowed_depth:
    """Find depth across windows in one sample."""
    input:
        bams = collect_runs_from_sample,  # From variant_calling.smk
        bams_idx = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_sample(wildcards))),
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
    conda: "../envs/delly.yaml"
    shell: """
        mosdepth {params.prefix} <(samtools merge {input.bam} -o -) \
            --by {params.window_size} \
            -n \
            --fast-mode \
            -t {threads} \
        """

rule merge_windowed_depth:
    """Merge windowed depth into one file for viewing joint figure."""
    input:
        #bed = expand(config["results"] + "coverage/mosdepth/{sample}.regions.bed.gz", sample=SAMPLES),
        bed = expand(config["results"] + "coverage/mosdepth/{sample}.regions.bed.gz",
            sample=SAMPLES),
    output:
        #merged_bed = config["results"] + "coverage/mosdepth/{name}.merged_with_dup.bed",
        merged_bed = config["results"] + "coverage/mosdepth/{name}.merged.bed",
    shell: """
        for BED in {input.bed}; do \
            SAMPLE=$(echo $BED | cut -d "." -f 1 | rev | cut -d "/" -f 1 | rev); \
            gunzip -c $BED \
            | awk -v SAMPLE=$SAMPLE 'BEGIN {{OFS="\\t"}} {{print SAMPLE, $0}}' \
            >> {output.merged_bed}; \
        done \
        """

## for figure --> notebooks/coverage.ipynb
