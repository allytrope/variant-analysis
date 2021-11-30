"""Contain rules for determining relatedness among individuals.
These rules will generally be used after those in `variant_calling.smk`.
"""

rule estimate_relatedness_lcmlkin:
    """Estimate relatedness between individuals using maximum likelihood."""
    input: vcf=config["lcmlkin"]["vcf"],
           founders=config["lcmlkin"]["founders"]
    output: config["output_path"] + "joint-call/Mmul_8.recalibrated_joint-call2.relate"
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
    output: ped=config["output_path"] + "plink/{dataset}.ped",
            map=config["output_path"] + "plink/{dataset}.map"
    shell: "vcftools \
                --vcf {input.vcf} \
                --plink \
                --out {wildcards.dataset}"

'''
rule create_binary_plink_files:
    """Create .bed and associated PLINK files."""
    input: "resources/"
    output: expand("{output_path}plink/data.{ext}", output_path=config["output_path"], ext=["bed", "bim", "fam"])
    shell: "plink \
                --vcf {input} #or maybe --file \
                --make-bed"
'''

rule recode_plink_files:
    input: ped=config["output_path"] + "plink/{dataset}.ped",
           map=config["output_path"] + "plink/{dataset}.map"
    output: ped=config["output_path"] + "plink/{dataset}_recode12.ped",
            map=config["output_path"] + "plink/{dataset}_recode12.map"
    wildcard_constraints: dataset="\w+"
    shell: "plink \
                --file  " + config["output_path"] + " plink/{wildcards.dataset} \
                --out " + config["output_path"] + " plink/{wildcards.dataset}_recode12 \
                --recode12"

rule admixture:
    """Estimate relatedness."""
    input: ped=config["output_path"] + "plink/{dataset}_recode12.ped",
           map=config["output_path"] + "plink/{dataset}_recode12.map"
    output: q=config["output_path"] + "relatedness/{dataset}.2.Q",
            p=config["output_path"] + "relatedness/{dataset}.2.P"
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
    output: config["output_path"] + "plink/{dataset}.pop"
    shell: "ls"
'''

rule supervised_admixture:
    """Estimate relatedness while taking into account known ancestry."""
    input: ped=config["output_path"] + "plink/{dataset}_recode12.ped",
           map=config["output_path"] + "plink/{dataset}_recode12.map",
           population=config["output_path"] + "plink/{dataset}.pop" # MAKE .POP FILE
    output: q=config["output_path"] + "relatedness/{dataset}.2.Q",
            p=config["output_path"] + "relatedness/{dataset}.2.P"
    params: ancestral_populations=2
    threads: 24
    shell: "admixture {input.ped} {params.ancestral_populations} \
                --supervised {input.population} \
                -j{threads}; \
            mv {wildcards.dataset}_recode12.{params.ancestral_populations}.Q {output.q}; \
            mv {wildcards.dataset}_recode12.{params.ancestral_populations}.P {output.p}"

# Under development
rule split_by_chrom:
    """Split data by chromosome."""
    input: ped=config["output_path"] + "plink/{dataset}.ped",
           map=config["output_path"] + "plink/{dataset}.map"
    #output: expand("{output_path}plink/{dataset}.{chr}.ped", output_path=config["output_path"], chr=list(range(1, 23)) + ['X']),
    #        expand("{output_path}plink/{dataset}.{chr}.map", output_path=config["output_path"], chr=list(range(1, 23)) + ['X'])
    output: config["output_path"] + "plink/{dataset}.{chr}.bed",
            config["output_path"] + "plink/{dataset}.{chr}.bim",
            config["output_path"] + "plink/{dataset}.{chr}.fam"
    shell: "plink \
                --file " + config["output_path"] + " plink/{wildcards.dataset}.{chr} \
                --chrom {chrom} \
                --out " + config["output_path"] + " plink/{wildcards.dataset}.{chr} \
                --make-bed"

rule haplotype_estimation:
    """Estimate haplotypes. Implementation assumes a constant recombination rate between SNPs."""
    input: bed=config["output_path"] + "plink/{dataset}.{chr}.bed",
           bim=config["output_path"] + "plink/{dataset}.{chr}.bim",
           fam=config["output_path"] + "plink/{dataset}.{chr}.fam"
    output: haps=config["output_path"] + "haplotypes/{dataset}.{chr}.shapeit.haps",
            sample=config["output_path"] + "haplotypes/{dataset}.{chr}.shapeit.sample"
    log: config["output_path"] + "haplotypes/{dataset}.{chr}.log"
    params: effective_size = config["population"]["effective_size"]
    threads: 24
    shell: "shapeit \
                --input-bed {input.bed} {input.bim} {input.fam} \
                --output-max {output.haps} {output.sample} \
                --output-log {log} \
                --efective-size {params.effective_size} \
                --thread {threads}"

rule prephasing:
    """Prep for imputation."""
    input: config["output_path"] + "haplotypes/{dataset}.{chr}"
    output: config["output_path"] + "haplotypes/{dataset}.{chr}.prephasing"
    params: map="",
            gens="",
            boundaries="1 5e6",  # IMPUTE2 suggests < 10 Mb
            effective_size = config["population"]["effective_size"]
    shell: "impute2 \
                -prephase_g \
                -m {params.map} \
                -g {params.gens} \
                -int {params.boundaries} \
                -Ne {params.effective_size} \
                -o {output}"

rule create_legend:
    """Create .legend file"""
    pass

rule impute:
    """Genotype imputation."""
    input: map=config["output_path"] + "haplotypes/{dataset}/{chr}.map",
           haps=config["output_path"] + "haplotypes/{dataset}/{chr}.haps",
           legend=config["output_path"] + "haplotypes/{dataset}/{chr}.legend",
           prephasing=config["output_path"] + "haplotypes/{dataset}/{chr}.prephasing"
           strand=""
    output: config["output_path"] + "haplotypes/{dataset}/{chr}"
    params: k_hap=500  #ADD VALUE TO DETERMINE NUMBER OF REFERENCE HAPLOTYPES
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

rule aggregate_haplotyped:
    """Combine haplotypes."""
    input: expand("{output_path}haplotypes/{dataset}.{chr}.phased.haps", output_path=config["output_path"], chr=list(range(1, 23)) + ['X']),
           expand("{output_path}haplotypes/{dataset}.{chr}.phased.sample", output_path=config["output_path"], chr=list(range(1, 23)) + ['X']),
    output: config["output_path"] + "aggregated"  # CHANGE
    shell: ""
