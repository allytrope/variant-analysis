"""Contain workflow for performing variant analysis from paired-end FASTQ files."""
configfile: "config.yaml"

# Path and sample names for input .fastq files
import os
FULL_SAMPLE_PATHS = [os.path.split(path) for path in config["data"]]
SAMPLE_PATHS, SAMPLE_NAMES = zip(*FULL_SAMPLE_PATHS)

rule all:
    input: expand("{output_path}bam/{sample}.bam", output_path=config["output_path"], sample=SAMPLE_NAMES)

# Prepare reference genome
rule download_ref:
    """Download reference genome."""
    output: config["ref"]
    shell: "rsync rsync://" + config["ref_address"] + " {output}"

rule index_ref:
    """Create BWA index files for reference genome."""
    input: config["ref"]
    output: multiext(config["ref"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    shell: "bwa index {input}"

def sample_path(sample):
    """Find path for corresponding sample name."""
    idx = SAMPLE_NAMES.index(sample)
    return SAMPLE_PATHS[idx] + "/"

# Alignment
rule align:
    """Align .fastq sequence to reference genome."""
    input: ref=config["ref"],
           idxs=multiext(config["ref"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
           read_1=lambda wildcards: sample_path("{sample}".format(sample=wildcards.sample)) + "{sample}.R1.fastq.gz",
           read_2=lambda wildcards: sample_path("{sample}".format(sample=wildcards.sample)) + "{sample}.R2.fastq.gz"
    output: config["output_path"] + "bam/{sample}.bam"
    shell: "bwa mem {input.ref} {input.read_1} {input.read_2} | \
            samtools view -Sb - > {output}"

# Post-alignment cleanup
rule remove_duplicates:
    """Mark duplicates with decimal value 1024."""
    input: config["output_path"] + "bam/{sample}.bam"
    output: mdup=config["output_path"] + "no_dupl_bam/{sample}.bam",
            metrics=config["output_path"] + "no_dupl_bam/{sample}.metrics.txt"
    conda: "envs/gatk.yaml"
    shell: "gatk MarkDuplicateSpark \
                -I {input} \
                -O {output.mdup} \
                -M {output.metrics} \
                --remove-duplicates true"

rule base_recalibration:
    """Create recalibration table for correcting systemic error in base quality scores."""
    input: ref=config["ref"],
           bam=config["output_path"] + "no_dupl_bam/{sample}.bam",
           known_variants=config["vcf"]
    output: config["output_path"] + "BQSR/{sample}.BQSR.recal"
    conda: "envs/gatk.yaml"
    shell: "gatk BaseRecalibrator \
                -R {input.ref} \
                -I {input.bam} \
                --known-sites {input.known_variants} \
                -O {output}"

rule apply_base_recalibration:
    """Correct systemic error in base quality scores."""
    input: ref=config["ref"],
           bam=config["output_path"] + "no_dupl_bam/{sample}.bam",
           recal=config["output_path"] + "BQSR/{sample}.BQSR.recal"
    output: config["output_path"] + "cleaned_bam/{sample}.bam"
    conda: "envs/gatk.yaml"
    shell: "gatk ApplyBQSR \
                -R {input.ref} \
                -I {input.bam} \
                --bqsr-recal-file {input.recal} \
                -O {output}"

# Variant calling
rule call_variants:
    """Call variants to make VCF file."""
    input: ref=config["ref"],
           bam=config["output_path"] + "cleaned_bam/{sample}.bam"
    output: config["output_path"] + "vcf/{sample}.g.vcf"
    conda: "envs/gatk.yaml"
    shell: 'gatk --java-options "-Xmx4g" HaplotypeCaller \
                -R {input.ref} \
                -I {input.bam} \
                -O {output} \
                -ERC GVCF'

rule create_sample_map:
    """Create sample map that contains names and paths to all VCFs to be used in consolidate rule."""
    input: expand("{output_path}vcf/{sample}.g.vcf", output_path=config["output_path"], sample=SAMPLE_NAMES)
    output: config["output_path"] + "/cohort.sample_map"
    run:
        import os
        with open("cohort.sample_map", "w") as out:
            for path, sample in FULL_SAMPLE_PATHS:
                string = sample + "\t" + config["output_path"] + "vcf/" + sample + ".vcf.gz" + "\n"
                out.write(string)


# Consolidation of GVCFs
rule consolidate:
    """Combines GVCFs into GenomicsDB (a datastore)."""
    input: expand("{output_path}{sample}.g.vcf", output_path=config["output_path"], sample=SAMPLE_NAMES), config["output_path"] + "cohort.sample_map"
    output: config["output_path"] + "db/genomicsdb"
    conda: "envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx4g -Xms4g' GenomicsDBImport \
                --sample-name-map cohort.sample_map \
                --genomicsdb-workspace-path db \
                --genomicsdb-update-workspace-path db"

rule joint_call_cohort:
    """Use GenomicsDB to jointly call a VCF file."""
    input: ref=config["ref"],
           db=config["output_path"] + "db"
    output: config["output_path"] + "joint-call.vcf.gz"
    conda: "envs/gatk.yaml"
    shell: "gatk --java-options '-Xmx4g' GenotypeGVCFs \
                -R {input.ref} \
                -V gendb://{input.db} \
                -O {output}"

# Variant quality score recalibration
rule calls_recalibration:
    """Build a recalibration model."""
    input: ref=config["ref"],
           vcf=config["output_path"] + "joint-call.vcf.gz",
           truth=config["truth_vcf"],
           training=config["training_vcf"]
    output: recal=config["output_path"] + "VQSR/VQSR.recal",
            tranches=config["output_path"] + "VQSR/output.tranches"
    conda: "envs/gatk.yaml"
    shell: "gatk VariantRecalibrator \
                -R {input.ref} \
                -V {input.vcf} \
                -O {output.recal} \
                --resource source1,known=false,truth=true,training=true,prior=10.0:{input.truth} \
                --resource source2,known=false,truth=false,training=true,prior=10.0:{input.training} \
                #-an add anotations from VCF \
                --tranches-file {output.tranches}"

rule apply_variant_recalibration:
    """Apply recalibration model."""
    input: ref=config["ref"],
           vcf=config["output_path"] + "joint-call.vcf.gz",
           recal=config["output_path"] + "VQSR/VQSR.recal",
           tranches=config["output_path"] + "VQSR/output.tranches"
    output: config["output_path"] + "recalibrated_joint-call.vcf.gz"
    conda: "envs/gatk.yaml"
    shell: "gatk ApplyVQSR \
                -R {input.ref} \
                -V {input.vcf} \
                -O {output} \
                --recal-file {input.recal} \
                --tranches-file {input.tranches}"
