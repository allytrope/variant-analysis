"""Contain rules for determining relatedness among individuals.
These rules will generally be used after those in `variant_calling.smk`.
"""

rule estimate_relatedness_lcmlkin:
    """Estimate relatedness between individuals using maximum likelihood."""
    input: vcf=config["lcmlkin"]["vcf"],
           founders=config["lcmlkin"]["founders"]
    output: config["output_path"] + "joint-call/Mmul_8.recalibrated_joint-call2.relate"
    threads: 8 # max of 8
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
'''
# Under development
rule split_by_chrom:
    """Split data by chromosome."""
    input: ped=config["output_path"] + "plink/{dataset}.ped",
           map=config["output_path"] + "plink/{dataset}.map"
    output: expand("{output_path}plink/{dataset}.{chr}.ped", output_path=config["output_path"], chr=list(range(1, 23)) + ['X']),
            expand("{output_path}plink/{dataset}.{chr}.map", output_path=config["output_path"], chr=list(range(1, 23)) + ['X'])
    shell: "for i in {{{{1..20}},X}}; do \
                echo '-B {dataset}.$1 \
                    -O'"
'''

'''for chr in {1..22}; do \
plink \
--file {} \
--chr {} \
--recode \
--out {}; \
done"'''

rule haplotype_estimation:
    """Estimate haplotypes. Implementation assumes a constant recombination rate between SNPs."""
    input: ped=config["output_path"] + "plink/{dataset}.{chr}.ped",
           map=config["output_path"] + "plink/{dataset}.{chr}.map"
    output: haps=config["output_path"] + "haplotypes/{dataset}.{chr}.phased.haps",
            sample=config["output_path"] + "haplotypes/{dataset}.{chr}.phased.sample"
    log: config["output_path"] + "haplotypes/{dataset}.{chr}.log"
    params: effective_size = config["shapeit"]["effective_population_size"]
    threads: 24
    shell: "shapeit \
                --input-ped {input.ped} {input.map}\
                --output-max {output.haps} {output.sample} \
                --output-log {log} \
                --efective-size {params.effective_size} \
                --thread {threads}"

'''
rule impute:
    """Genotype imputation."""
    input: ""
    output: ""
    shell: "impute2"
'''
