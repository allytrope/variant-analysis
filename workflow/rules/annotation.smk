"""Annotate VCFs."""


rule bcsq:
    """Add BCSQ tag for consequences."""
    input:
        bcf = config["results"] + "haplotypes/SHAPEIT5/{dataset}.{subset}.{mode}.chr{chr}.bcf",
        csi = config["results"] + "haplotypes/SHAPEIT5/{dataset}.{subset}.{mode}.chr{chr}.bcf.csi",
        gff = config["resources"] + "annotations/Macaca_mulatta.Mmul_10.110.gff3.gz",        #gff = config["gtf3"], 
        ref_fasta = config["ref_fasta"],
        #gff_idx = config["gtf3"] + ".tbi",
    output:
        bcf = config["results"] + "genotypes/annotated/bcsq/{dataset}.{subset}.{mode}.chr{chr}.bcf",
    conda:  "../envs/common.yaml",
    shell: """
        bcftools csq {input.bcf} \
            -f {input.ref_fasta} \
            -g {input.gff} \
            -Ob \
            -o {output.bcf} \
        """

rule count_bcsq_consequences:
    """Parse BCSQ to TSV. This keep only sample for each individual and counts up the type of consequences."""
    input:
        bcf = config["results"] + "genotypes/annotated/bcsq/{dataset}.{subset}.{mode}.chr{chr}.bcf",
        samples = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.dedup_samples.list",
    output:
        tsv = config["results"] + "genotypes/annotated/bcsq/{dataset}.{subset}.{mode}.chr{chr}.tsv",
    shell: """
        bcftools view {input.bcf} \
            -S {input.samples} \
            -Ou \
        | bcftools +split-vep \
            -f '%CHROM:%POS\t%gene\t%biotype\t%transcript\t%Consequence\t%INFO/AC\t%INFO/AN\n' \
            -s 'worst' \
        | grep -v '@' \
        | cut -f 5 \
        | sort \
        | uniq -c \
        > {output.tsv}
        """

rule early_vep:
    """Run VEP at an early step than usual."""
    input:
        vcf = config["results"] + "joint_call/polyallelic/{dataset}.chr{chr}.vcf.gz",
        fasta = config["ref_fasta"],
        gtf = config["gtf3"],
        gtf_idx = config["gtf3"] + ".tbi",
    output:
        annot = config["results"] + "joint_call/annotated/{dataset}.chr{chr}.vcf.gz",
    params:
        #species = "macaca_mulatta",
        species = "callithrix_jacchus",
    threads: 1
    resources: nodes = 1
    conda:  "../envs/ensembl.yaml",
    #             --compress_output bgzip \
    shell: """
        bcftools view {input.vcf} \
            -Ov \
        | vep \
            --fasta {input.fasta} \
            --format vcf \
            --gtf {input.gtf} \
            --output_file STDOUT \
            --species {params.species} \
            --vcf \
        | bcftools view \
            -Oz \
            -o {output.annot} \
        """

rule vep_SNV:
    """Annotate using ENSEMBL VEP."""
    input:
        bcf = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.bcf",
        fasta = config["ref_fasta"],
        gtf = config["gtf3"],
        gtf_idx = config["gtf3"] + ".tbi",
    output:
        annot = config["results"] + "genotypes/annotated/{dataset}.{subset}.{mode}.chr{chr}.vcf",
        stats = config["results"] + "genotypes/annotated/{dataset}.{subset}.{mode}.chr{chr}.stats.html",
    params:
        species = "macaca_mulatta",
        output_format = "--vcf"
    threads: 1
    resources: nodes = 1
    conda:  "../envs/ensembl.yaml",
    #             --compress_output bgzip \
    # | bcftools view \
    # -Ob \
    # -o {output.annot} \
    shell: """
        bcftools view {input.bcf} \
            -Ov \
        | vep \
            --fasta {input.fasta} \
            --format vcf \
            --gtf {input.gtf} \
            --output_file STDOUT \
            --stats_file {output.stats} \
            --species {params.species} \
            {params.output_format} \
            -o {output.annot} \
        """
use rule vep_SNV as vep_SNV_to_TSV with:
    output:
        annot = config["results"] + "genotypes/annotated/{dataset}.{subset}.{mode}.chr{chr}.tsv",
        stats = config["results"] + "genotypes/annotated/{dataset}.{subset}.{mode}.chr{chr}.stats.html",
    params:
        species = "macaca_mulatta",
        output_format = ""  # This defaults VEP to output to TSV format

# Should no longer be needed since I'm inheriting the above rule
# rule vep_SNV_to_TSV:
#     """Annotate using ENSEMBL VEP to TSV."""
#     input:
#         bcf = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.bcf",
#         fasta = config["ref_fasta"],
#         gtf = config["gtf3"],
#         gtf_idx = config["gtf3"] + ".tbi",
#     output:
#         annot = config["results"] + "genotypes/annotated/{dataset}.{subset}.{mode}.chr{chr}.tsv",
#     params:
#         species = "macaca_mulatta",
#     threads: 1
#     resources: nodes = 1
#     conda:  "../envs/ensembl.yaml",
#     #             --compress_output bgzip \
#     shell: """
#         bcftools view {input.bcf} \
#             -Ov \
#         | vep \
#             --fasta {input.fasta} \
#             --format vcf \
#             --gtf {input.gtf} \
#             --output_file {output.annot} \
#             --species {params.species} \
#         """

# Condensed rules into one
rule annotate_and_filter:
    """Annotate and then filter to only high or moderate impact SNVs."""
    input:
        bcf = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.bcf",
        fasta = config["ref_fasta"],
        gtf = config["gtf3"],
        gtf_idx = config["gtf3"] + ".tbi",
    output:
        annot = config["results"] + "genotypes/impactful/{dataset}.{subset}.{mode}.chr{chr}.bcf",
        stats = config["results"] + "genotypes/impactful/{dataset}.{subset}.{mode}.chr{chr}.html",
    params:
        species = "macaca_mulatta",
    threads: 1
    resources: nodes = 1
    conda:  "../envs/ensembl.yaml",
    #             --compress_output bgzip \
    shell: """
        bcftools view {input.bcf} \
            -Ov \
        | vep \
            --fasta {input.fasta} \
            --format vcf \
            --gtf {input.gtf} \
            --output_file STDOUT \
            --species {params.species} \
            --stats_file {output.stats} \
            --vcf \
        | bcftools +split-vep \
            -c 'IMPACT' \
        | bcftools view -i 'IMPACT="HIGH" | IMPACT="MODERATE"' \
            -Ob \
            -o {output.annot} \
        """

# rule filter_high_and_moderate_impact:
#     input:
#         annot = config["results"] + "genotypes/annotated/{dataset}.{subset}.{mode}.chr{chr}.bcf",
#     output:
#         annot = config["results"] + "genotypes/impactful/{dataset}.{subset}.{mode}.chr{chr}.bcf",
#     shell: """
#         bcftools +split-vep {input.annot} \
#             -c 'IMPACT' \
#         | bcftools view -i 'IMPACT="HIGH" | IMPACT="MODERATE"' \
#             -Ob \
#             -o {output.annot} \
#         """


# rule create_liftover_ref_dict:
#     """Create .dict file for compressed reference genome."""
#     #wildcard_constraints: fasta = "fasta|fna|fa"
#     input: fasta = "/master/abagwell/variant-analysis/resources/rhesus/annotations/Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
#     output: dict = "/master/abagwell/variant-analysis/resources/rhesus/annotations/Homo_sapiens.GRCh38.dna.toplevel.dict",
#     threads: 1
#     resources: nodes = 1
#     conda: "../envs/gatk.yaml"
#     shell: "gatk CreateSequenceDictionary \
#             -R {input.fasta}"

rule vep_clinvar:
    """Annotate liftover using ClinVar, which requires a human genome."""
    input:
        bcf = config["results"] + "liftover/{liftover}/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz",
        clinvar = config["resources"] + "annotations/clinvar_20260208.vcf.gz",
        clinvar_idx = config["resources"] + "annotations/clinvar_20260208.vcf.gz.tbi",
        fasta = config["ref_fasta"],
        #gtf = config["gtf3"],
        #gtf_idx = config["gtf3"] + ".tbi",
    output:
        annot = config["results"] + "liftover/{liftover}/annotated/VEP/{dataset}.{subset}.{mode}.chr{chr}.tsv",
    threads: 1
    resources: nodes = 1
    conda:  "../envs/ensembl.yaml",
    # --gtf {input.gtf} \
    shell: """
        bcftools view {input.bcf} \
            -Ov \
        | vep \
            --fasta {input.fasta} \
            --format vcf \
            --output_file {output.annot} \
            --species homo_sapiens \
            --assembly GRCh38 \
            --custom file={input.clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN \
        """


# rule split_BND_breakpoints:
#     """Create new lines to separate a BND into its breakpoints."""
#     input:
#         vcf = "/master/abagwell/variant-analysis/results/rhesus/structural_variants/SVs/merged/{dataset}.genotyped.pass.vcf.gz",
#     output:
#         vcf = "/master/abagwell/variant-analysis/results/rhesus/structural_variants/SVs/merged/{dataset}.genotyped.pass.split_BNDs.vcf.gz",
#     # Note that this still gives this message, but otherwise works:
#     # [W::vcf_parse_format] Extreme FORMAT/RC value encountered and set to missing at 6:20788029
#     run:
#         ## This code WOULD work, but variant.POS does not have a setter (although variant.POS does)
#         from cyvcf2 import VCF, Writer

#         # Open I/O
#         vcf = VCF("U42_WGS_WES.genotyped.pass.vcf.gz")
#         w = Writer("U42_WGS_WES.genotyped.pass.split_BNDs.vcf.gz", vcf)

#         # Iterate through variants
#         for variant in vcf:
#             w.write_record(variant)
#             if variant.INFO.get("CHR2"): # FIX
#                 chr2 = variant.INFO.get("CHR2")
#                 pos2 = variant.INFO.get("POS2")
#                 end = variant.INFO.get("END")

#                 variant.CHROM = chr2
#                 #variant.POS = pos2
#                 variant.set_pos(pos2 - 1)


#                 replaced = str(variant).replace(str(end), str(pos2 + 1)).strip()

#                 replaced_variant = w.variant_from_string(replaced)


#                 # Write new record
#                 w.write_record(replaced_variant)


#         # Close I/O
#         w.close(); vcf.close()

rule split_BND_breakpoints:
    """Create new lines to separate a BND into its breakpoints."""
    input:
        vcf = "/master/abagwell/variant-analysis/results/rhesus/structural_variants/delly/merged/{dataset}.genotyped.pass.vcf.gz",
    output:
        vcf = "/master/abagwell/variant-analysis/results/rhesus/structural_variants/delly/merged/{dataset}.genotyped.pass.split_BNDs.vcf.gz",
    # Note that this still gives this message, but otherwise works:
    # [W::vcf_parse_format] Extreme FORMAT/RC value encountered and set to missing at 6:20788029
    resources: nodes = 1
    conda: "../envs/delly2.yaml"
    script: "workflow/scripts/split_BND_breakpoints.py"


# run:
#     with open(input.vcf, "r") as f, open(output.vcf, "w") as g:
#         for line in f:
#             if line.startswith("#"):
#                 g.write(line)
#             else:
#                 g.write(line)
#                 columns = line.split("\t")
#                 if "CHR2" in columns[7]:
#                     chr2 = [field for field in columns[7].split(";") if "CHR2" in field][0].split("=")[1]
#                     pos2 = [field for field in columns[7].split(";") if "POS2" in field][0].split("=")[1]
#                     columns[0] = chr2
#                     columns[1] = pos2
#                     new_row = "\t".join(columns)
#                     g.write(new_row)


rule vep_SV:
    """Annotate using ENSEMBL VEP."""
    input:
        vcf = "/master/abagwell/variant-analysis/results/rhesus/structural_variants/SVs/merged/{dataset}.genotyped.pass.split_BNDs.vcf.gz",
        fasta = config["ref_fasta"],
        gtf = config["gtf3"],
        gtf_idx = config["gtf3"] + ".tbi",
    output:
        annot = "/master/abagwell/variant-analysis/results/rhesus/structural_variants/annotated/{dataset}.genotyped.pass.split_BNDs.vcf.gz",
    params:
        species = "macaca_mulatta",
    threads: 1
    resources: nodes = 1
    conda:  "../envs/ensembl.yaml",
    shell: """
        vep \
            --input_file {input.vcf} \
            --fasta {input.fasta} \
            --species {params.species} \
            --gtf {input.gtf} \
            --vcf \
            --compress_output bgzip \
            --output_file {output.annot} \
        """

# rule vep_SNV:
#     """Annotate using ENSEMBL VEP."""
#     input:
#         #vcf = "/master/abagwell/variant-analysis/results/rhesus/haplotypes/SHAPEIT5_merged/{dataset}.SNP.chr{chr}.bcf",
#         #bcf = config["results"] + "gwas/SHAPEIT5_WGS/bcfs/{dataset}.SNP.chr{chr}.bcf",
#         vcf = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.vcf.gz",
#         tbi = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
#         fasta = config["ref_fasta"],
#         gtf = config["gtf3"],
#         gtf_idx = config["gtf3"] + ".tbi",
#     output:
#         annot = config["results"] + "genotypes/annotated/{dataset}.SNP.chr{chr}.vcf.gz",
#     params:
#         species = "macaca_mulatta",
#     threads: 1
#     resources: nodes = 1
#     conda:  "../envs/ensembl.yaml",
#     shell: """
#         vep \
#             --input_file {input.vcf} \
#             --fasta {input.fasta} \
#             --species {params.species} \
#             --gtf {input.gtf} \
#             --vcf \
#             --compress_output bgzip \
#             --output_file {output.annot} \
#         """

# Working on
rule igv:
    """Plot variants."""
    input:
    output:
    shell: """
        create_report test/data/variants/variants.bed \
            --genome hg38 \
            --flanking 1000 \
            --info-columns GENE TISSUE TUMOR COSMIC_ID GENE SOMATIC \
            --tracks test/data/variants/variants.bed test/data/variants/recalibrated.bam \
            --output examples/example_genome.html
        """


# rule snp_summary:
#     """Gives a summary of SNP data from .vcf file."""
#     input: vcf=config["results"] + "joint_call/recalibrated_joint_call.vcf.gz"
#     output: config["results"] + "joint_call/snp_summary.txt"
#     shell: "vep \
#     --custom /home/flow/ensembl-vep/Mmul_8.0.1.92.chr.gff3.gz,,gff \
#     --fasta Mmul_8.0.1.chromosome.fa \
#     --gff /home/flow/ensembl-vep/Mmul_8.0.1.92.chr.gff3.gz \
#     --input_file {input.vcf} \
#     --output_file {output} \
#     --stats_text "


# Download cache for species of interest first
# rule snp_summary:
#     """Gives a summary of SNP data from .vcf file."""
#     input: vcf = config["results"] + "joint_call/{workspace}_recalibrated.vcf.gz",
#            ref = config["results"] + "ref/ref_genome.fna.gz",
#     output: txt = config["results"] + "joint_call/{workspace}_snp_summary.txt",
#             html = config["results"] + "joint_call/{workspace}_snp_summary.html",
#             warnings = config["results"] + "joint_call/{workspace}_snp_summary_warnings.txt",
#     params: species = config["ref"]["species"],
#     threads: 1
#     conda: "../envs/common.yaml"
#     shell: "vep \
#     --input_file {input.vcf} \
#     --fasta {input.ref} \
#     --output_file {output.txt} \
#     --stats_text {output.html} \
#     --warning_file {output.warnings} \
#     --species {params.species} \
#     --cache"

rule bcf_stats:
    """Compute summary stats."""
    input:
        bcfs = expand(config["results"] + "genotypes/pass/{{dataset}}.{{subset}}.SNP.chr{chr}.bcf",
            chr=AUTOSOMES),
    output:
        stats = config["results"] + "genotypes/annotated/{dataset}.{subset}.SNP.autosomal.stats.tsv",
    threads: 2
    resources: nodes = 1
    conda:  "../envs/common.yaml",
    shell: """
        bcftools stats <(bcftools concat {input.bcfs}) > {output.stats} \
        """