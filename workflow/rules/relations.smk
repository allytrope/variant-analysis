"""Contain rules for determining relatedness among individuals.
These rules will generally be used after those in `variant_calling.smk`.
"""

rule estimate_relatedness_lcmlkin:
    """Estimate relatedness between individuals using maximum likelihood."""
    input: vcf=config["lcmlkin"]["vcf"],
           founders=config["lcmlkin"]["founders"]
    output: config["results"] + "joint-call/recalibrated_joint-call.relate"
    threads: 8  # max of 8
    shell: "lcmlkin \
                -i {input.vcf} \
                -o {output} \
                -g all \
                -l phred \
                -u {input.founders} \
                -t {threads}"

rule create_plink_files:
    """Create PLINK .ped and .map files."""
    input: vcf=config["lcmlkin"]["vcf"]
    output: ped=config["results"] + "plink/{dataset}.ped",
            map=config["results"] + "plink/{dataset}.map"
    shell: "vcftools \
                --vcf {input.vcf} \
                --plink \
                --out {wildcards.dataset}"

'''
rule create_binary_plink_files:
    """Create .bed and associated PLINK files."""
    input: "resources/"
    output: expand("{results}plink/data.{ext}", results=config["results"], ext=["bed", "bim", "fam"])
    shell: "plink \
                --vcf {input} #or maybe --file \
                --make-bed"
'''

rule recode_plink_files:
    """Convert base readings from {A,C,T,G} to {1,2}."""
    input: ped=config["results"] + "plink/{dataset}.ped",
           map=config["results"] + "plink/{dataset}.map"
    output: ped=config["results"] + "plink/recode12/{dataset}_recode12.ped",
            map=config["results"] + "plink/recode12/{dataset}_recode12.map"
    shell: "plink \
                --file  " + config["results"] + " plink/{wildcards.dataset} \
                --out " + config["results"] + " plink/recode12/{wildcards.dataset}_recode12 \
                --recode12"

rule admixture:
    """Estimate relatedness."""
    input: ped=config["results"] + "plink/recode12/{dataset}_recode12.ped",
           map=config["results"] + "plink/recode12/{dataset}_recode12.map"
    output: q=config["results"] + "relatedness/{dataset}.2.Q",
            p=config["results"] + "relatedness/{dataset}.2.P"
    params: ancestral_populations=2
    threads: 24
    shell: "admixture {input.ped} {params.ancestral_populations} \
                -j{threads}; \
            mv {wildcards.dataset}_recode12.{params.ancestral_populations}.Q {output.q}; \
            mv {wildcards.dataset}_recode12.{params.ancestral_populations}.P {output.p}"

# Under development
'''
rule create_pop:
    """Create admixture .pop file."""
    output: config["results"] + "plink/{dataset}.pop"
    shell: "ls"
'''

rule supervised_admixture:
    """Estimate relatedness while taking into account known ancestry."""
    input: ped=config["results"] + "plink/recode12/{dataset}_recode12.ped",
           map=config["results"] + "plink/recode12/{dataset}_recode12.map",
           population=config["results"] + "plink/{dataset}.pop" # MAKE .POP FILE
    output: q=config["results"] + "relatedness/{dataset}.2.Q",
            p=config["results"] + "relatedness/{dataset}.2.P"
    params: ancestral_populations=2
    threads: 24
    shell: "admixture {input.ped} {params.ancestral_populations} \
                --supervised {input.population} \
                -j{threads}; \
            mv {wildcards.dataset}_recode12.{params.ancestral_populations}.Q {output.q}; \
            mv {wildcards.dataset}_recode12.{params.ancestral_populations}.P {output.p}"

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

rule create_gen:
    """Create .gen file for IMPUTE."""
    input: ped=config["results"] + "plink/split/{dataset}.chr{chr}.ped",
           map=config["results"] + "plink/split/{dataset}.chr{chr}.map"
    output: config["results"] + "haplotypes/split/{dataset}.chr{chr}.gen"
    shell: "gtool -P --ped {input.ped} --map {input.map} --og {output}"

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
