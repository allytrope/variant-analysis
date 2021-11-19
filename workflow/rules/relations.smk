"""Contain rules for determining relatedness among individuals.
These rules will generally take place after those in `variant_calling.smk`.
"""

'''
rule vcf2geno:
    input:
    output:
    shell:
'''

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
'''
rule create_ped:
    """Create PLINK .ped file."""
    input: script=ancient("scripts/create_ped.awk"),
           pedigree=config["plink"]["pedigree"]
    output: config["output_path"] + "plink/GbS_2020.ped"
    shell: "awk -f {input.script} {input.pedigree} > {output}"

rule create_map:
    """Create PLINK .map file."""
    input: script=ancient("scripts/create_map.awk"),
           vcf=config["lcmlkin"]["vcf"]
    output: config["output_path"] + "plink/GbS_2020.map"
    shell: "egrep -v '^#' {input.vcf} | awk -f {input.script} > {output}"
'''

rule create_ped_files:
    """Create PLINK .ped and .map files."""
    input: vcf=config["lcmlkin"]["vcf"]
    output: ped=config["output_path"] + "plink/{dataset}.ped",
            map=config["output_path"] + "plink/{dataset}.map"
    shell: "vcftools --vcf {input.vcf} --plink --out {wildcards.dataset}"

'''
rule create_bed:
    """Create .bed and associated PLINK files."""
    input: "resources/"
    output: expand("{output_path}plink/data.{ext}", output_path=config["output_path"], ext=["bed", "bim", "fam"])
    shell: "plink \
            --vcf {input} #or maybe --file \
            --make-bed"
'''

rule recode_ped_files:
    input: ped=config["output_path"] + "plink/{name}.ped",
           map=config["output_path"] + "plink/{name}.map"
    output: ped=config["output_path"] + "plink/{name}_recode12.ped",
            map=config["output_path"] + "plink/{name}_recode12.map"
    wildcard_constraints: name="\w+"
    shell: "plink --file  " + config["output_path"] + " plink/{wildcards.name} --recode12 --out " + config["output_path"] + " plink/{wildcards.name}_recode12"

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
            mv {wildcards.dataset}.{params.ancestral_populations}.Q {output.q}; \
            mv {wildcards.dataset}.{params.ancestral_populations}.P {output.p}"

'''
rule create_genotypefile:
    """Create .geno file for LASER."""
    # from .pef or .vcf?
    input: config["lcmlkin"]["vcf"]
    output: "output.geno", "output.site"
    shell: "vcf2geno \
            --inVcf \
            --out output"

rule create_sitefile:
    """Crate .site file for LASER."""
    input:
    output:
    shell:

rule laser:
    """Estimate individual ancestry background."""
    input: ""
    output: config["output_path"] + "relations/laser.coord"
    shell: ""

rule estimate_relatedness_seekin:
    """Estimate relatedness with SEEKIN, which is better for closely related individuals and admixture."""
    input: ""
    output: ""
    #threads: 
    shell: "seekin kinship \
            -i {input} \
            -c {input}"
'''

'''
rule shapeit4:
    """Estimate haplotypes."""
    input: vcf="", map=config["shapeit4"]["map"]
    output: config["output_path2"] + "joint-call/Mmul_8.phased.vcf"
    log: config["output_path2"] + "joint-call/phased.log"
    threads: 24
    shell: "shapeit4 \
            --input {input.vcf} \
            --map {input.map} \
            --region rrrrr \
            --output {output} \
            --log {log} \
            --sequencing \
            --thread {threads}"
'''