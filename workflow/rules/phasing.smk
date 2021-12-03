"""Contain rules for phasing. The PLINK .ped and .map inputs come from `relations.smk`."""

rule split_by_chrom:
    """Split data by chromosome."""
    input: ped=config["results"] + "plink/{dataset}.ped",
           map=config["results"] + "plink/{dataset}.map"
    output: config["results"] + "plink/split/{dataset}.chr{chr}.ped",
            config["results"] + "plink/split/{dataset}.chr{chr}.map"
    shell: "plink \
                --file " + config["results"] + "plink/{wildcards.dataset} \
                --chr {wildcards.chr} \
                --out " + config["results"] + "plink/split/{wildcards.dataset}.chr{wildcards.chr} \
                --recode \
                --noweb"

rule create_gen:
    """Create .gen file for IMPUTE."""
    input: ped=config["results"] + "plink/split/{dataset}.chr{chr}.ped",
           map=config["results"] + "plink/split/{dataset}.chr{chr}.map"
    output: config["results"] + "haplotypes/split/{dataset}.chr{chr}.gen"
    log: config["log"] + "gtool/gen/{dataset}.chr{chr}.log"
    shell: "gtool -P \
           --ped {input.ped} \
           --map {input.map} \
           --og {output} \
           --log {log}"

rule haplotype_estimation:
    """Estimate haplotypes. Implementation assumes a constant recombination rate between SNPs."""
    input: ped=config["results"] + "plink/split/{dataset}.chr{chr}.ped",
           map=config["results"] + "plink/split/{dataset}.chr{chr}.map"
    output: haps=config["results"] + "haplotypes/split/{dataset}.chr{chr}.shapeit.haps",
            sample=config["results"] + "haplotypes/split/{dataset}.chr{chr}.shapeit.sample"
    params: effective_size = config["population"]["effective_size"]
    threads: 24
    shell: "shapeit \
                --input-ped {input.ped} {input.map} \
                --output-max {output.haps} {output.sample} \
                --efective-size {params.effective_size} \
                --thread {threads}"

rule prephasing:
    """Prep for imputation."""
    input: map=config["results"] + "haplotypes/split/{dataset}.chr{chr}.map",
           gen=config["results"] + "haplotyptes/split/{dataset}.chr{chr}.gen"
    output: config["results"] + "haplotypes/split/{dataset}.chr{chr}.prephasing.impute2"
    params: boundaries="1 5e6",  # IMPUTE2 suggests < 10 Mb
            effective_size = config["population"]["effective_size"]
    shell: "impute2 \
                -prephase_g \
                -m {input.map} \
                -g {input.gen} \
                -int {params.boundaries} \
                -Ne {params.effective_size} \
                -o {output}"

#rule create_legend:
#    """Create .legend file"""
#    pass

rule impute:
    """Genotype imputation."""
    input: map=config["results"] + "haplotypes/split/{dataset}.chr{chr}.map",
           haps=config["results"] + "haplotypes/split/{dataset}.chr{chr}.haps",
           legend=config["results"] + "haplotypes/split/{dataset}.chr{chr}.legend",
           prephasing=config["results"] + "haplotypes/split/{dataset}.chr{chr}.prephasing.impute2",
           strand=""
    output: config["results"] + "haplotypes/split/{dataset}.chr{chr}.phased.impute2"
    params: k_hap=500,  #ADD VALUE TO DETERMINE NUMBER OF REFERENCE HAPLOTYPES
            effective_size=config["population"]["effective_size"],  # IMPUTE2 suggests 20000
            boundaries="1 5e6"  # IMPUTE2 suggests < 10 Mb
    shell: "impute2 \
                -use_prephased_g \
                -m {input.map} \
                -h {input.haps} \
                -l {input.legend} \
                -known_haps_g {input.prephasing} \
                -strand_g {input.strand} \
                -k_hap {params.k_hap} \
                -Ne {params.effective_size} \
                -int {params.boundaries} \
                -o {output} \
                -phase"
'''
rule aggregate_haplotyped:
    """Combine haplotypes."""
    input: expand("{results}haplotypes/{dataset}.chr{chr}.phased.haps", results=config["results"], dataset=?, chr=list(range(1, 23)) + ['X']),
           expand("{results}haplotypes/{dataset}.chr{chr}.phased.sample", results=config["results"], dataset=?, chr=list(range(1, 23)) + ['X'])
    output: config["results"] + "aggregated"  # CHANGE
    shell: "ls"
'''
