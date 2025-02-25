"""Rules related to variant calling of SNVs."""

## Variant calling

rule call_variants:
    """Call variants from a single chromosome to make VCF file."""
    input:
        ref = config["ref_fasta"],
        fai = config["ref_fasta"] + ".fai",
        dict = ".".join(config["ref_fasta"].split(".")[:-2]) + ".dict",
        #bam = config["results"] + "alignments/merged/{sample}.bam",
        bams = lambda wildcards: [f'{config["results"]}alignments/recalibrated/{seqsample_library_run}.bam' for seqsample_library_run in collect_runs_from_library(wildcards, full=True)],
        bais = lambda wildcards: [f'{config["results"]}alignments/recalibrated/{seqsample_library_run}.bam.bai' for seqsample_library_run in collect_runs_from_library(wildcards, full=True)],
        #bam_i = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_library(wildcards, full=True))),
    params:
        bams = lambda wildcards, input: list(map(lambda bam: "-I " + bam, input.bams)),
    output:
        vcf = config["results"] + "gvcf/{chr}/{seq}{indiv_id}_{library}.chr{chr}.g.vcf.gz",
        #tbi = config["results"] + "gvcf/{chr}/{seq}{indiv_id}_{library}.chr{chr}.g.vcf.gz.tbi",
    log:
        config["results"] + "gvcf/{chr}/{seq}{indiv_id}_{library}.chr{chr}.log",
    conda: "../envs/gatk.yaml"
    threads: 2  # 4 is default. Also see https://hpc.nih.gov/training/gatk_tutorial/haplotype-caller.html for recommended threads.
    resources: nodes = 2
    shell: """
        gatk --java-options '-Xmx16g' HaplotypeCaller \
            -L {wildcards.chr} \
            -R {input.ref} \
            {params.bams} \
            -O {output.vcf} \
            -ERC GVCF \
            --native-pair-hmm-threads {threads} \
            2> {log} \
        """

rule call_variants_merge_chromosomes:
    """Merge called variants from single chromosomes to make VCF file."""
    wildcard_constraints:
        indiv_id = r"[A-Za-z0-9]+",
        seq = r"WGS|WES|GBS|AMP",
        run = r"[A-Za-z0-9_-]+",
    input:
        vcfs = lambda wildcards: expand(config["results"] + "gvcf/{chr}/{seq}{indiv_id}_{library}.chr{chr}.g.vcf.gz",
            chr=CHROMOSOMES,
            seq = wildcards.seq,
            indiv_id = wildcards.indiv_id,
            library = wildcards.library),
    output:
        vcf = config["results"] + "gvcf/{seq}{indiv_id}_{library}.g.vcf.gz",
        #tbi = config["results"] + "gvcf/{seq}{indiv_id}_{library}.g.vcf.gz.tbi",
    conda: "../envs/bio.yaml"
    threads: 1
    resources: nodes = 1
    shell: """
        bcftools concat {input.vcfs} \
            -o {output.vcf} \
            -Oz \
        """

# Consolidation of GVCFs
rule create_sample_map:
    """Create sample map that contains names and paths to all VCFs to be used in consolidate rule.
    Note: This output file will need to be deleted if changing what will be added in the consolidate rule."""
    input:
        gvcfs = lambda wildcards: expand("{results}gvcf/{chr}/{sample}.chr{chr}.g.vcf.gz",
            results=config["results"],
            sample=SAMPLES,
            chr=wildcards.chr),
    output:
        #sample_map = temp(config["results"] + "db/{dataset}.sample_map"),
        sample_map = config["results"] + "db/sample-maps/{dataset}.chr{chr}.sample-map",
    threads: 1
    resources: nodes = 1
    run:
        import json
        with open(f'/master/abagwell/variant-analysis/results/rhesus/db/{wildcards.dataset}/1/callset.json') as f:
            data = json.load(f)
            samples_already_in_datastore = [sample['sample_name'] for sample in data['callsets']]

        with open(output.sample_map, "w") as sample_map:
            for gvcf in input.gvcfs:
                sample = gvcf.split("/")[-1].split(".")[0]
                if sample not in samples_already_in_datastore:
                    sample_map.write(f"{sample}\t{gvcf}\n")

rule consolidate:
    """Combine the chromosomes of .g.vcf files into GenomicsDB datastore."""
    input:
        contigs = config["resources"] + "ref_fna/chromosomes.list",
        gvcfs = lambda wildcards: expand("{results}gvcf/{chr}/{sample}.chr{chr}.g.vcf.gz{ext}",
        #gvcfs = lambda wildcards: expand("{results}gvcf/{sample}.g.vcf.gz{ext}",
            ext=["", ".tbi"],
            results=config["results"],
            sample=SAMPLES,
            chr=wildcards.chr),
        sample_map = config["results"] + "db/sample-maps/{dataset}.chr{chr}.sample-map",  # gVCFs are referenced in this file
    output:
        touch(config["results"] + "db/{dataset}/completed/{chr}.txt"),
    params:
        db = config["results"] + "db/{dataset}/{chr}",
        # Higher value requires more memory and number of file descriptor able to be used at the same time.
        # Attempted with 6, but ran out of memory. Once even 4 was too much. Though did work with 4 when I had fewer samples.
        parallel_intervals = 3,
    threads: 2  # Just for opening multiple .vcf files at once.
    resources: nodes = 2
    conda: "../envs/gatk.yaml"
    # CONTIGS=$(awk 'BEGIN {{ORS = ","}} {{print $0}}' {input.contigs}); \
    shell: """
        if [ -d {params.db} ]; \
        then WORKSPACE_FLAG="genomicsdb-update-workspace-path"; \
        else WORKSPACE_FLAG="genomicsdb-workspace-path"; \
        fi; \
        gatk --java-options '-Xmx8g' GenomicsDBImport \
            --$WORKSPACE_FLAG {params.db} \
            --intervals {wildcards.chr} \
            --sample-name-map {input.sample_map} \
            --batch-size 50 \
            --genomicsdb-shared-posixfs-optimizations true \
            --reader-threads {threads} \
            --max-num-intervals-to-import-in-parallel {params.parallel_intervals} \
        """

# rule joint_call_trio:
#     """Use GenomicsDB to jointly call a VCF file."""
#     input:
#         ref = config["ref_fasta"],
#         # Note: The actual text file isn't what is required, but the datastore directory.
#         # This .txt file is, however, is created only after the datastore has finished being built.
#         db = config["results"] + "db/{dataset}/completed/{contig}.txt",
#     output:
#         vcf = config["results"] + "joint_call/polyallelic/{dataset}.chr{contig}.vcf.gz",
#         #tbi = config["results"] + "joint_call/polyallelic/{dataset}.chr{chr}.vcf.gz.tbi",
#     params:
#         db = config["results"] + "db/{dataset}/{contig}",
#     threads: 1
#     resources: nodes = 1
#     conda: "../envs/gatk.yaml"
#     shell: """
#         gatk --java-options '-Xmx16g' GenotypeGVCFs \
#             -R {input.ref} \
#             -V gendb://{params.db} \
#             -O {output.vcf} \
#             -L {wildcards.contig} \
#         """

## Jointly call variants
rule join_call_cohort:
    """Use GenomicsDB to jointly call a VCF file."""
    input:
        ref = config["ref_fasta"],
        # Note: The actual text file isn't what is required, but the datastore directory.
        # This .txt file is, however, is created only after the datastore has finished being built.
        db = config["results"] + "db/{dataset}/completed/{contig}.txt",
    output:
        vcf = config["results"] + "joint_call/polyallelic/{dataset}.chr{contig}.vcf.gz",
        #tbi = config["results"] + "joint_call/polyallelic/{dataset}.chr{chr}.vcf.gz.tbi",
    params:
        db = config["results"] + "db/{dataset}/{contig}",
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    shell: """
        gatk --java-options '-Xmx16g' GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{params.db} \
            -O {output.vcf} \
            -L {wildcards.contig} \
        """

rule joint_call_cohort_all_sites:
    """Use GenomicsDB to jointly call a VCF file. This version call all invariant sites as well."""
    input:
        ref = config["ref_fasta"],
        # Note: The actual text file isn't what is required, but the datastore directory.
        # This .txt file is, however, is created only after the datastore has finished being built.
        db = config["results"] + "db/{dataset}/completed/{contig}.txt",
    output:
        vcf = config["results"] + "joint_call/all_sites/{dataset}.chr{contig}.vcf.gz",
        #tbi = config["results"] + "joint_call/polyallelic/{dataset}.chr{chr}.vcf.gz.tbi",
    params:
        db = config["results"] + "db/{dataset}/{contig}",
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    shell: """
        gatk --java-options '-Xmx16g' GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{params.db} \
            -O {output.vcf} \
            -L {wildcards.contig} \
            -all-sites \
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
        #-e'type{params.equality}"snp"' \
        #equality = lambda wildcards: "=" if wildcards.mode == "indel" else "!=",
        equality_option = lambda wildcards: """-e'type="snp"'""" if wildcards.mode == "indel" else (
            """-e'type!="snp"'""" if wildcards.mode == "SNP" else ""
            ),

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
            {params.equality_option} \
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
