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
                --make-bed --noweb"
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
