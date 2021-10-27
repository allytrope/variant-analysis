"""Contain all rules for performing variant analysis from paired-end FASTQ files."""
configfile: "config.yaml"

# Prepare reference genome
rule download_ref:
    """Download reference genome."""
    output: config["ref"]
    shell: "rsync rsync://" + config["ref_address"]
rule index_ref:
    """Create BWA index files for reference genome."""
    input: config["ref"]
    output: multiext(config["ref"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    shell: "bwa index {input}"

# Alignment
rule align:
    """Align .fastq sequence to reference genome."""
    input: ref=config["ref"],
           idxs=multiext(config["ref"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
           read_1=config["datasets_path"] + "/{dataset}_1.fq.gz",
           read_2=config["datasets_path"] + "/{dataset}_2.fq.gz"
    output: "bam/{dataset}.bam"
    shell: "bwa mem {input.ref} {input.read_1} {input.read_2} \
            | samtools view -Sb - > {output}"

# Post-alignment cleanup
rule mark_duplicates:
    """Mark duplicates with decimal value 1024."""
    input: "bam/{dataset}.bam"
    output: mdup="marked_duplicates/{dataset}.bam",
            metrics="marked_duplicates/{dataset}.metrics.txt"
    container: "docker://broadinstitute/gatk"
    shell: "gatk MarkDuplicateSpark \
                -I {input} \
                -O {output.mdup} \
                -M {output.metrics}"

# Variant calling
rule call_variants:
    """Call variants to make VCF file."""
    input: ref=config["ref"], bam="marked_duplicates/{dataset}.bam"
    output: "vcf/{dataset}.vcf"
    container: "docker://broadinstitute/gatk"
    shell: 'gatk --java-options "-Xmx4g" HaplotypeCaller \
                -R {input.ref} \
                -I {input.bam} \
                -O {output} \
                -ERC GVCF'

rule create_sample_map:
    """Create sample map that contains names and paths to all VCFs to be used in consolidate rule."""
    input: ["vcf/{dataset}.vcf".format(dataset=dataset) for dataset in config["datasets"]]
    params: datasets=[dataset for dataset in shell("find fastq -type f -name '*.fq'", iterable=True)]
    output: "cohort.sample_map"
    run:
        import os
        with open("cohort.sample_map", "w") as out:
            for location in params.datasets:
                base_file = os.path.basename(location)
                base_name = os.path.splitext(base_file)[0]
                string = base_name + "\t" + "fastq/" + base_name + ".vcf" + "\n"
                out.write(string)

# Consolidation of GVCFs
rule consolidate:
    """Combines GVCFs into GenomicsDB (a datastore)."""
    input: ["vcf/{dataset}.vcf".format(dataset=dataset) for dataset in config["datasets"]], "cohort.sample_map"
    output: "db/genomicsdb"
    container: "docker://broadinstitute/gatk"
    shell: "gatk --java-options '-Xmx4g -Xms4g' GenomicsDBImport \
                --sample-name-map cohort.sample_map \
                --genomicsdb-workspace-path db \
                --genomicsdb-update-workspace-path db"

rule joint_call_cohort:
    """Use GenomicsDB to jointly call a VCF file."""
    input: ref=config["ref"], db="db"
    output: "joint-call.vcf.gz"
    container: "docker://broadinstitute/gatk"
    shell: "gatk --java-options '-Xmx4g' GenotypeGVCFs \
                -R {input.ref} \
                -V gendb://{input.db} \
                -O {output}"

# Variant quality score recalibration
rule calls_recalibration:
    """Build a recalibration model."""
    input: ref=config["ref"], vcf="joint-call.vcf.gz"
    output: "output.recal"
    container: "docker://broadinstitute/gatk"
    shell: "gatk VariantRecalibrator \
                -R {input.ref} \
                -V {input.vcf} \
                -O {output.recal}"

rule apply_recalibration:
    """Apply recalibration model."""
    input: ref=config["ref"], vcf="joint-call.vcf.gz"
    output: "recalibrated_joint-call.vcf.gz"
    container: "docker://broadinstitute/gatk"
    shell: "gatk ApplyVQSR \
                -R {input.ref} \
                -V {input.vcf} \
                -O {output}"
