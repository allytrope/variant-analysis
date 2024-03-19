"""Rules related to variant calling of SNVs."""

## Variant calling

rule call_variants:
    """Call variants to make VCF file."""
    input:
        ref = config["ref_fasta"],
        fai = config["ref_fasta"] + ".fai",
        dict = ".".join(config["ref_fasta"].split(".")[:-2]) + ".dict",
        #bam = config["results"] + "alignments/merged/{sample}.bam",
        bams = collect_runs_from_sample,
        bams_idx = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_sample(wildcards)))
        #bam_idx = config["results"] + "alignments/merged/{sample}.bam.bai",
    params:
        bams = lambda wildcards, input: list(map(lambda bam: "-I " + bam, input.bams)),
    output:
        vcf = config["results"] + "gvcf/{sample}.g.vcf.gz",
        tbi = config["results"] + "gvcf/{sample}.g.vcf.gz.tbi",
    conda: "../envs/gatk.yaml"
    threads: 9  # 4 is default
    resources: nodes = 9
    # -I {input.bam} \
    shell: """
        gatk --java-options '-Xmx16g' HaplotypeCaller \
            -R {input.ref} \
            {params.bams} \
            -O {output.vcf} \
            -ERC GVCF \
            --native-pair-hmm-threads {threads} \
        """

# Consolidation of GVCFs
rule create_sample_map:
    """Create sample map that contains names and paths to all VCFs to be used in consolidate rule.
    Note: This output file will need to be deleted if changing what will be added in the consolidate rule."""
    input:
        gvcfs = expand("{results}gvcf/{sample}.g.vcf.gz",
            results=config["results"],
            sample=SAMPLES),
    output:
        #sample_map = temp(config["results"] + "db/{dataset}.sample_map"),
        sample_map = config["results"] + "db/{dataset}.sample-map",
    threads: 1
    resources: nodes = 1
    run:
        with open(output.sample_map, "w") as sample_map:
            for gvcf in input.gvcfs:
                sample = gvcf.split("/")[-1].split(".")[0]
                sample_map.write(f"{sample}\t{gvcf}\n")

rule list_chromosomes:
    """Find all chromosomes from headers of reference genome.
    
     Each line of the output file is just one chromosome's name.
     The search only keeps numbered chromosomes (not those prefixed with "chr" or any other letters)
     as well as X, Y, and MT. Unplaced contigs are ignored."""
    input:
        ref = config["ref_fasta"],
    output:
        config["results"] + "db/chromosomes.list",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        zcat {input.ref} | grep "^>" | cut -c 2- | cut -d " " -f 1 | grep -x -E "^[0-9]+|X|Y|MT" > {output}
        """

rule consolidate:
    """Combine the chromosomes of .g.vcf files into GenomicsDB datastore."""
    input:
        contigs = config["results"] + "db/chromosomes.list",
        sample_map = config["results"] + "db/{dataset}.sample-map",
    output:
        config["results"] + "db/created_{dataset}.txt",
    params:
        db = config["results"] + "db/{dataset}",
        # Higher value requires more memory and number of file descriptor able to be used at the same time.
        # Attempted with 6, but ran out of memory. Once even 4 was too much. Though did work with 4 when I had fewer samples.
        parallel_intervals = 3,  
    threads: 2  # Just for opening multiple .vcf files at once.
    resources: nodes = 2
    conda: "../envs/gatk.yaml"
    shell: """
        CONTIGS=$(awk 'BEGIN {{ORS = ","}} {{print $0}}' {input.contigs}); \
        if [ -d {params.db} ]; \
        then WORKSPACE_FLAG="genomicsdb-update-workspace-path"; \
        else WORKSPACE_FLAG="genomicsdb-workspace-path"; \
        fi; \
        gatk --java-options '-Xmx8g' GenomicsDBImport \
            --$WORKSPACE_FLAG {params.db} \
            --intervals {input.contigs} \
            --sample-name-map {input.sample_map} \
            --batch-size 50 \
            --genomicsdb-shared-posixfs-optimizations true \
            --reader-threads {threads} \
            --max-num-intervals-to-import-in-parallel {params.parallel_intervals} \
        && touch {output}
        """

## Jointly call variants
rule joint_call_cohort:
    """Use GenomicsDB to jointly call a VCF file."""
    input:
        ref = config["ref_fasta"],
        # Note: The actual text file isn't what is required, but the datastore directory.
        # This .txt file is, however, is created only after the datastore has finished being built.
        db = config["results"] + "db/created_{dataset}.txt",
    output:
        vcf = config["results"] + "joint_call/polyallelic/{dataset}.chr{chr}.vcf.gz",
        #tbi = config["results"] + "joint_call/polyallelic/{dataset}.chr{chr}.vcf.gz.tbi",
    params:
        db = config["results"] + "db/{dataset}",
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    shell: """
        gatk --java-options '-Xmx16g' GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{params.db} \
            -O {output.vcf} \
            -L {wildcards.chr} \
        """

## Split SNPs and indels
rule biallelics_by_mode:
    """Split into SNP- or indel-only .vcf. Then keeps only biallelic sites."""
    input:
        vcf = config["results"] + "joint_call/polyallelic/{dataset}.chr{chr}.vcf.gz",
        tbi = config["results"] + "joint_call/polyallelic/{dataset}.chr{chr}.vcf.gz.tbi",
        ref_fasta = config["ref_fasta"],
    output:
        split = config["results"] + "joint_call/biallelic/{dataset}.{mode}.chr{chr}.vcf.gz",
    params:
        equality = lambda wildcards: "=" if wildcards.mode == "indel" else "!=",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    # 1) Separate multiallelics into different lines
    # 2) Take only SNPs or indels
    # 3) Merge multiallelics back into same lines
    # 4) Keep only biallelics
    # Alternative:
    # bcftools view {input.vcf} \
    # -M2 \
    # -v snps \
    # -Oz \
    # -o {output.split} \
    shell: """
        bcftools norm {input.vcf} \
            -m-any \
            --fasta-ref {input.ref_fasta} \
            -Ou \
        | bcftools view \
            -e'type{params.equality}"snp"' \
            -Ou \
        | bcftools norm \
            -m+any \
            -Ou \
        | bcftools view \
            -M2 \
            -m2 \
            -Oz \
            -o {output.split} \
        """
