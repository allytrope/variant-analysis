"""Contain rules for imputation. Runs after phasing.smk."""

## Using SHAPEIT

rule make_scaffold:
    """Create scaffold for SHAPEIT4."""
    input:
        vcf = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
        fam = config["results"] + "haplotypes/pedigree/all_samples.trios_only.tsv",
    output:
        scaffold = config["results"] + "haplotypes/scaffolds/{dataset}.{mode}.chr{chr}.vcf.gz",
    threads: 1
    resources: nodes = 1
    shell: """
        makeScaffold \
            --gen {input.vcf} \
            --fam {input.fam} \
            --reg {wildcards.chr} \
            --out {output.scaffold} \
        """

# rule shapeit4_imputation_ref:
#     """Haplotype estimation and imputation for the reference VCF, that is to say, WGS samples only."""
#     input:
#         vcf = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz",
#         csi = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz.csi",
#         #scaffold = config["results"] + "haplotypes/scaffolds/{dataset}.{mode}.chr{chr}.vcf.gz",
#         #scaffold_csi = config["results"] + "haplotypes/scaffolds/{dataset}.{mode}.chr{chr}.vcf.gz.csi",
#     output:
#         phased = config["results"] + "haplotypes/SHAPEIT4/{dataset}.{mode}.chr{chr}.vcf.gz",
#     # When using SHAPIT4.1 or greater, --map not required (though surely still helpful).
#     # chr appears to still be mandatory even if there is only one chromosome in file.
#     params:
#         PS = 0.0002,  # 0.0001 is recommended value by SHAPEIT4
#     log: config["results"] + "haplotypes/SHAPEIT4/log/12_Indian_12_Chinese.chr{chr}.log",
#     threads: 8
#     resources: nodes = 8
#     conda: "../envs/shapeit4.yaml"
#     shell: """
#         shapeit4 \
#             --input {input.vcf} \
#             --region {wildcards.chr} \
#             --sequencing \
#             --use-PS {params.PS} \
#             --output {output.phased} \
#             --log {log} \
#             --thread {threads} \
#         """

# rule shapeit4_imputation_query:
#     """Haplotype estimation and imputation."""
#     input:
#         vcf = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz",
#         csi = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz.csi",
#         #map = config["resources"] + "genetic_map/chr{chr}.cM.genetic_map",
#         scaffold = config["results"] + "haplotypes/scaffolds/{dataset}.{mode}.chr{chr}.vcf.gz",
#         scaffold_csi = config["results"] + "haplotypes/scaffolds/{dataset}.{mode}.chr{chr}.vcf.gz.csi",
#     output:
#         phased = config["results"] + "haplotypes/SHAPEIT4_{seq}/{dataset}.{mode}.chr{chr}.vcf.gz",
#     # When using SHAPIT4.1 or greater, --map not required (though surely still helpful).
#     # chr appears to still be mandatory even if there is only one chromosome in file.
#     params:
#         PS = 0.0002,  # 0.0001 is recommended value by SHAPEIT4
#     log: config["results"] + "haplotypes/SHAPEIT4_new/log/{dataset}.{mode}.chr{chr}.log",
#     threads: 8
#     resources: nodes = 8
#     conda: "../envs/shapeit4.yaml"
#     shell: """
#         shapeit4 \
#             --input {input.vcf} \
#             --scaffold {input.scaffold} \
#             --region {wildcards.chr} \
#             --sequencing \
#             --use-PS {params.PS} \
#             --output {output.phased} \
#             --log {log} \
#             --thread {threads} \
#         """

# rule shapeit4_imputation_ref_vcf:
#     """Haplotype estimation and imputation for the reference VCF used for determining admixture."""
#     input:
#         #vcf = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz",
#         vcf = "/master/abagwell/variant-analysis/resources/rhesus/ref_vcf/12_Indian_12_Chinese.vcf.gz",
#         #csi = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz.csi",
#         csi = "/master/abagwell/variant-analysis/resources/rhesus/ref_vcf/12_Indian_12_Chinese.vcf.gz.csi",
#         #scaffold = config["results"] + "haplotypes/scaffolds/{dataset}.{mode}.chr{chr}.vcf.gz",
#         #scaffold_csi = config["results"] + "haplotypes/scaffolds/{dataset}.{mode}.chr{chr}.vcf.gz.csi",
#     output:
#         #phased = config["results"] + "haplotypes/SHAPEIT4/with_scaffold/{dataset}.{mode}.chr{chr}.vcf.gz",
#         phased = "/master/abagwell/variant-analysis/resources/rhesus/ref_vcf/12_Indian_12_Chinese.phased.chr{chr}.vcf.gz",
#     # When using SHAPIT4.1 or greater, --map not required (though surely still helpful).
#     # chr appears to still be mandatory even if there is only one chromosome in file.
#     params:
#         PS = 0.0001,  # Recommended value by SHAPEIT4
#     log: config["results"] + "haplotypes/SHAPEIT4/log/12_Indian_12_Chinese.chr{chr}.log",
#     threads: 8
#     resources: nodes = 8
#     conda: "../envs/shapeit4.yaml"
#     shell: """
#         shapeit4 \
#             --input {input.vcf} \
#             --region {wildcards.chr} \
#             --sequencing \
#             --use-PS {params.PS} \
#             --output {output.phased} \
#             --log {log} \
#             --thread {threads} \
#         """

rule add_annotations:
    """Adds FORMAT annotations that were removed during processing with SHAPEIT4."""
    input: 
        vcf = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz",
        phased = config["results"] + "haplotypes/SHAPEIT4_{seq}/{dataset}.{mode}.chr{chr}.vcf.gz",
    output:
        annotated = config["results"] + "haplotypes/SHAPEIT4_{seq}/annotated/{dataset}.{mode}.chr{chr}.vcf.gz",
    threads: 1 
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools annotate {input.vcf} \
            -a {input.phased} \
            -c FORMAT/GT \
            -o {output.annotated} \
            -Oz \
        """

rule shapeit5_pedigree:
    """Generate pedigree file for SHAPEIT5. Three tab-delimited columns.
    Only includes duos and trios. And fills unknown parent with NA."""
    input:
        tsv = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.tsv",
    output:
        fam = config["results"] + "haplotypes/pedigree/{dataset}.duos_and_trios.fam",
    threads: 1
    resources: nodes = 1
    conda: "../envs/shapeit5.yaml"
    shell: """
        sed 's/\\t\\t/\\tNA\\t/g;s/\\t$/\\tNA/g' {input.tsv} \
        | grep -v NA$'\\t'NA \
        > {output.fam} \
        """

rule shapeit5_imputation_ref:
    """Haplotype estimation and imputation of reference VCF, that is to say, WGS samples."""
    input:
        vcf = config["results"] + "genotypes/pass/WGS/{dataset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "genotypes/pass/WGS/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
        fam = config["results"] + "haplotypes/pedigree/{dataset}.duos_and_trios.fam",
        #map = config["resources"] + "genetic_map/chr{chr}.cM.genetic_map",
    output:
        phased = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.{mode}.chr{chr}.bcf",
    log:
        config["results"] + "haplotypes/SHAPEIT5_WGS/log/{dataset}.{mode}.chr{chr}.log",
    threads: 8
    resources: nodes = 8
    conda: "../envs/shapeit5.yaml"
    shell: """
        SHAPEIT5_phase_common \
            --input {input.vcf} \
            --pedigree {input.fam} \
            --region {wildcards.chr} \
            --output {output.phased} \
            --log {log} \
            --thread {threads} \
        """
#            --map {input.map} \

rule shapeit5_imputation:
    """Haplotype estimation and imputation of WES samples using the WGS imputation as the reference VCF."""
    input:
        vcf = config["results"] + "genotypes/pass/WES/{dataset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "genotypes/pass/WES/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
        fam = config["results"] + "haplotypes/pedigree/{dataset}.duos_and_trios.fam",
        ref = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.{mode}.chr{chr}.bcf",
        #map = config["resources"] + "genetic_map/chr{chr}.cM.genetic_map",
    output:
        phased = config["results"] + "haplotypes/SHAPEIT5_WES/{dataset}.{mode}.chr{chr}.bcf",
    log:
        config["results"] + "haplotypes/SHAPEIT5_WES/log/{dataset}.{mode}.chr{chr}.log",
    threads: 8
    resources: nodes = 8
    conda: "../envs/shapeit5.yaml"
    shell: """
        SHAPEIT5_phase_common \
            --input {input.vcf} \
            --pedigree {input.fam} \
            --reference {input.ref} \
            --region {wildcards.chr} \
            --output {output.phased} \
            --log {log} \
            --thread {threads} \
        """
#            --map {input.map} \

## Using Beagle

rule beagle_imputation:
    """Impute using Beagle."""
    input:
        vcf = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz",
    output:
        vcf = config["results"] + "haplotypes/Beagle5/{dataset}.{mode}.chr{chr}.vcf.gz",
    # log:
    #     log = config["results"] + "haplotypes/Beagle5/{dataset}.{mode}.chr{chr}.vcf.gz.log",
    params:
        jar = "/master/abagwell/tools/beagle.22Jul22.46e.jar",
        out = lambda wildcards, output: ".".join(output.vcf.split(".")[:-2]),
    threads: 24
    resources: nodes = 24
    # Beagle5 not available in conda
    conda: "../envs/gatk.yaml"  # Using this enviornment because it already has Java8
    shell: """
        java -Xmx50g -jar {params.jar} \
            gt={input.vcf} \
            out={output.out} \
            nthreads={threads} \
        """


## Using GLIMPSE

# rule shapeit4_WGS_imputation:
#     """Haplotype estimation and imputation."""
#     input:
#         vcf = config["results"] + "haplotypes/whatshap/WGS/{dataset}.{mode}.chr{chr}.vcf.gz",
#         csi = config["results"] + "haplotypes/whatshap/WGS/{dataset}.{mode}.chr{chr}.vcf.gz.csi",
#         #map = config["resources"] + "genetic_map/chr{chr}.cM.genetic_map",
#         scaffold = config["results"] + "haplotypes/scaffolds/{dataset}.{mode}.chr{chr}.vcf.gz",
#         scaffold_csi = config["results"] + "haplotypes/scaffolds/{dataset}.{mode}.chr{chr}.vcf.gz.csi",
#     output:
#         phased = config["results"] + "haplotypes/SHAPEIT4_WGS/{dataset}.{mode}.chr{chr}.vcf.gz",
#     # When using SHAPIT4.1 or greater, --map not required (though surely still helpful).
#     # chr appears to still be mandatory even if there is only one chromosome in file.
#     params:
#         PS = 0.0001,  # 0.0001 is recommended value by SHAPEIT4
#     log: config["results"] + "haplotypes/SHAPEIT4_WGS/log/{dataset}.{mode}.chr{chr}.log",
#     threads: 8
#     resources: nodes = 8
#     conda: "../envs/shapeit4.yaml"
#     shell: """
#         shapeit4 \
#             --input {input.vcf} \
#             --scaffold {input.scaffold} \
#             --region {wildcards.chr} \
#             --sequencing \
#             --use-PS {params.PS} \
#             --output {output.phased} \
#             --log {log} \
#             --thread {threads} \
#         """

rule shapeit5_WGS_imputation:
    """Haplotype estimation and imputation."""
    input:
        vcf = config["results"] + "genotypes/pass/WGS/{dataset}.{mode}.chr{chr}.vcf.gz",
        csi = config["results"] + "genotypes/pass/WGS/{dataset}.{mode}.chr{chr}.vcf.gz.csi",
        #vcf = config["results"] + "genotypes/pass/WGS/{dataset}.{mode}.chr{chr}.vcf.gz",
        #csi = config["results"] + "genotypes/pass/WGS/{dataset}.{mode}.chr{chr}.vcf.gz.csi",
        fam = config["results"] + "haplotypes/pedigree/{dataset}.duos_and_trios.fam",
        #map = config["resources"] + "genetic_map/chr{chr}.cM.genetic_map",
        #scaffold = config["results"] + "haplotypes/scaffolds/{dataset}.{mode}.chr{chr}.vcf.gz",
        #scaffold_csi = config["results"] + "haplotypes/scaffolds/{dataset}.{mode}.chr{chr}.vcf.gz.csi",
    output:
        phased = config["results"] + "haplotypes/SHAPEIT5/WGS/{dataset}.{mode}.chr{chr}.bcf",
    log:
        config["results"] + "haplotypes/SHAPEIT5/WGS/log/{dataset}.{mode}.chr{chr}.log",
    threads: 8
    resources: nodes = 8
    conda: "../envs/shapeit5.yaml"
    shell: """
        SHAPEIT5_phase_common \
            --input {input.vcf} \
            --pedigree {input.fam} \
            --region {wildcards.chr} \
            --output {output.phased} \
            --log {log} \
            --thread {threads} \
        """

rule glimpse_chunk:
    """Chunk into regions for GLIMPSE imputation."""
    input:
        vcf = config["results"] + "haplotypes/whatshap/WES/{dataset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "haplotypes/whatshap/WES/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
        #map = config["results"] + "genetic_map/chr{chr}.cM.genetic_map",
    output:
        chunks = config["results"] + "haplotypes/glimpse/chunks/{dataset}.{mode}.chr{chr}.chunks.txt",
    log:
        config["results"] + "haplotypes/glimpse/chunks/{dataset}.{mode}.chr{chr}.chunks.log",
    threads: 24
    resources: nodes = 24
    conda: "../envs/glimpse.yaml"
    shell: """
        GLIMPSE2_chunk \
            --input {input.vcf} \
            --region {wildcards.chr} \
            --sequential \
            --output {output.chunks} \
            --log {log} \
            --threads {threads} \
        """
        # --map {input.map} \

# rule glimpse_reference:
#     """Format reference VCF for GLIMPSE use."""
#     input:
#         vcf = config["results"] + "haplotypes/SHAPEIT5/WGS/{dataset}.{mode}.chr{chr}.bcf",
#         csi = config["results"] + "haplotypes/SHAPEIT5/WGS/{dataset}.{mode}.chr{chr}.bcf.csi",
#         #map = expand(config["results"] + "genetic_map/chr{chr}.cM.genetic_map"),
#         chunks = config["results"] + "haplotypes/glimpse/chunks/{dataset}.{mode}.chr{chr}.chunks.txt",
#     output:
#         ref = config["results"] + "haplotypes/glimpse/ref/{dataset}.{mode}.chr{chr}.{region}.bin",
#     log:
#         log = config["results"] + "haplotypes/glimpse/ref/{dataset}.{mode}.chr{chr}.bin.log",
#     params:
#         out_prefix = ".".join(output.ref.split(".")[:-2]),
#     threads: 8
#     resources: nodes = 8
#     conda: "../envs/glimpse.yaml"
#     #             --map {input.map} \
#     shell: """
#         while IFS="" read -r LINE || [ -n "$LINE" ]; \
#         do \
#             IRG=$(echo $LINE | cut -d" " -f3);
#             ORG=$(echo $LINE | cut -d" " -f4);
#             GLIMPSE2_split_reference \
#                 --reference {input.vcf} \
#                 --input-region $IRG \
#                 --output-region $ORG \
#                 --output {output.out_prefix} \
#                 --threads {threads} \
#                 --log {log.log}; \
#         done < {input.chunks}
#         """

rule glimpse_reference_per_region:
    """Format reference VCF for GLIMPSE use."""
    input:
        vcf = config["results"] + "haplotypes/SHAPEIT5/WGS/{dataset}.{mode}.chr{chr}.bcf",
        csi = config["results"] + "haplotypes/SHAPEIT5/WGS/{dataset}.{mode}.chr{chr}.bcf.csi",
        #map = expand(config["results"] + "genetic_map/chr{chr}.cM.genetic_map"),
        chunks = config["results"] + "haplotypes/glimpse/chunks/{dataset}.{mode}.chr{chr}.chunks.txt",
    output:
        ref = config["results"] + "haplotypes/glimpse/ref/{dataset}.{mode}._{chr}_{start}_{end}.bin",
    log:
        log = config["results"] + "haplotypes/glimpse/ref/{dataset}.{mode}._{chr}_{start}_{end}.log",
    params:
        out_prefix = lambda wildcards, output: ".".join(output.ref.split(".")[:-2]) + ".",
    threads: 8
    resources: nodes = 8
    conda: "../envs/glimpse.yaml"
    #             --map {input.map} \
    shell: """
            ORG=$(grep {wildcards.chr}:{wildcards.start}-{wildcards.end} {input.chunks} | cut -f 4); \
            GLIMPSE2_split_reference \
                --reference {input.vcf} \
                --input-region {wildcards.chr}:{wildcards.start}-{wildcards.end} \
                --output-region $ORG \
                --output {params.out_prefix} \
                --threads {threads} \
                --log {log.log}; \
        """



# # Work in progress
# def collect_chr_intervals(wildcards, input):
#     """Collect all interval runs of GLIMPSE_split_reference."""
#     sample_runs = []
#     for run in SAMPLE_RUNS:
#         if wildcards.sample in run:
#             sample_runs.append(config["results"] + "alignments/recalibrated/" + run + ".bam")
#     return sample_runs
# rule collect_glimpse_reference:
#     """Call variants to make VCF file."""
#     input:
#         chunks = config["results"] + "haplotypes/glimpse/{dataset}.{mode}.chr{chr}.chunks.txt",
#         intervals = collect_chr_intervals,
#         bams_idx = lambda wildcards, : list(map(lambda bam: bam + ".bai", collect_chr_intervals(wildcards)))
#     params:
#         bams = lambda wildcards, input: list(map(lambda bam: "-I " + bam, input.bams)),
#     output:
#         vcf = 
#     conda: "../envs/gatk.yaml"
#     threads: 1
#     resources: nodes = 1
#     shell: """
#         """

# Will need to modify for when I use indiviuals with multiple BAMs
rule create_WES_bam_list:
    """Create list of BAM files for GLIMPSE2 phasing."""
    input:
        #config["results"] + "alignments/recalibrated/{sample_run}.bam",
        vcf = config["results"] + "haplotypes/whatshap/WES/{dataset}.{mode}.chr1.vcf.gz",
    output:
        bam_list = config["results"] + "haplotypes/glimpse/regions/{dataset}.{mode}.WES_bams.list",
    params:
        prefix = config["results"] + "alignments/recalibrated/",
        suffix = ".bam",
    shell: """
        for SAMPLE in $(bcftools query -l {input.vcf}); do \
            echo {params.prefix}${{SAMPLE}}{params.suffix}; \
        done > {output.bam_list} \
        """

rule glimpse_imputation:
    """Impute using GLIMPSE. For low-coverage."""
    input:
        bam_list = config["results"] + "haplotypes/glimpse/regions/{dataset}.{mode}.WES_bams.list",
        ref = config["results"] + "haplotypes/glimpse/ref/{dataset}.{mode}._{chr}_{start}_{end}.bin",
    output:
        vcf = config["results"] + "haplotypes/glimpse/regions/{dataset}.{mode}.chr{chr}_{start}_{end}.bcf",
    threads: 24
    resources: nodes = 24
    conda: "../envs/glimpse.yaml"
    shell: """
        GLIMPSE2_phase \
            --bam-list {input.bam_list} \
            --reference {input.ref} \
            --output {output.vcf} \
            --threads {threads} \
        """

def find_chunks(wildcards):
    """Return list of chunk file names for chromosome, based on the chunks.txt file."""
    chunk_files = []
    chunks = config["results"] + f"haplotypes/glimpse/chunks/{wildcards.dataset}.{wildcards.mode}.chr{wildcards.chr}.chunks.txt"
    with open(chunks, "r") as f:
        for line in f:
            pos = line.split("\t")[2]
            chrom = pos.split(":")[0]
            start = pos.split(":")[1].split("-")[0]
            end = pos.split(":")[1].split("-")[1]
            chunk_files.append(config["results"] + f"haplotypes/glimpse/regions/{wildcards.dataset}.{wildcards.mode}.chr{wildcards.chr}_{start}_{end}.bcf",)
    return chunk_files
rule glimpse_imputed_list:
    """Create list of imputed files."""
    input:
        chunks = config["results"] + "haplotypes/glimpse/chunks/{dataset}.{mode}.chr{chr}.chunks.txt",
        chunk_files = lambda wildcards: find_chunks(wildcards),
    output:
        txt = config["results"] + "haplotypes/glimpse/regions/{dataset}.{mode}.chr{chr}.file_names.txt",
    params:
        chunks = config["results"] + "haplotypes/glimpse/chunks/{dataset}.{mode}.chr{chr}.chunks.txt",
        prefix = config["results"] + "haplotypes/glimpse/regions/{dataset}.{mode}.chr",
        suffix = ".bcf",
    run:
        with open(output.txt, "w") as f:
            for line in input.chunk_files:
                f.write(line + "\n")

    
# """
# CHR=$(cut -f 2 {input.chunks})
# START=$(cut -f 3 {input.chunks} | cut -d ":" -f 2 | cut -d "-" -f 1)
# END=$(cut -f 3 {input.chunks} | cut -d ":" -f 2 | cut -d "-" -f 2)

# {params.prefix}${{CHR}}_${{START}}_${{END}}{params.suffix}

# """"

rule glimpse_ligate:
    """Combine phased regions into single chromosomes."""
    input:
        txt = config["results"] + "haplotypes/glimpse/regions/{dataset}.{mode}.chr{chr}.file_names.txt",
    output:
        vcf = config["results"] + "haplotypes/glimpse/chr/{dataset}.{mode}.chr{chr}.bcf",
    # log:
    #     log = config["results"] + "haplotypes/glimpse/chr/{dataset}.{mode}.chr{chr}.log",
    threads: 24
    resources: nodes = 24
    conda: "../envs/glimpse.yaml"
    shell: """
        GLIMPSE2_ligate \
            --input {input.txt} \
            --output {output.vcf} \
            --threads {threads} \
        """


## AlphaImpute2
# rule AlphaImpute2_impute:
#     """Impute using AlphaImpute2."""
#     input:
#         genotypes = 
#         pedigree =
#     output:
#         vcf = config["results"] + "haplotypes/glimpse/{dataset}.{mode}.chr{chr}.bcf",
#     params:
#         out_dir = lambda wildcards, output: ".".join(output.vcf.split(".")[:-1]),
#     threads: 16
#     resources: nodes = 16
#     conda: "../envs/alphaimpute2.yaml"
#     shell: """
#         AlphaImpute2 \
#             -genotypes {input.genotypes} \
#             -pedigree {input.pedigree} \
#             -out {params.out_dir} \
#             -maxthreads {threads} \
#         """