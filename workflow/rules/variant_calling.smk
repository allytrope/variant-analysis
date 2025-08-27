"""Rules related to variant calling of SNVs."""

## Variant calling

# rule variant_calling_bcftools:
#     """Variant calling with `bcftools`."""
#     shell: """
#         bcftools mpileup {input.bam} \
#             -Ou \
#             -f {input.ref_fasta} \
#         | bcftools call \
#             --ploidy 2 \
#             -m \
#             -v \
#             -o {output.bcf} \
#             -Ob \
#         """

rule variant_calling_freebayes:
    """Variant calling with `freebayes`."""
    shell: """
        freebayes {input.bams} \
            -f {input.ref_fasta} \
            -r {wildcards.chr} \
        """

rule call_variants:
    """Call variants from a single chromosome to make VCF file.
    Verify that rule is grouping by the appropriate field, using library."""
    input:
        ref = config["ref_fasta"],
        fai = config["ref_fasta"] + ".fai",
        dict = ".".join(config["ref_fasta"].split(".")[:-2]) + ".dict",
        bams = lambda wildcards: expand(config["results"] + "alignments/markdup/{sample}.bam",
            sample=collect_samples(
                fmt="{batch}/{seq}{indiv}_{library}_{flowcell_lane}",
                col="library",
                val=wildcards.library),
                # col="indiv",
                # val=wildcards.indiv),
        ),
        bais = lambda wildcards: expand(config["results"] + "alignments/markdup/{sample}.bam.bai",
            sample=collect_samples(
                fmt="{batch}/{seq}{indiv}_{library}_{flowcell_lane}",
                col="library",
                val=wildcards.library),
                # col="indiv",
                # val=wildcards.indiv),
        ),
    params:
        bams = lambda wildcards, input: list(map(lambda bam: "-I " + bam, input.bams)),
    output:
        vcf = config["results"] + "gvcf/{batch}/{seq}{indiv}_{library}.chr{chr}.g.vcf.gz",
        #tbi = config["results"] + "gvcf/{chr}/{seq}{indiv}_{library}.chr{chr}.g.vcf.gz.tbi",
    log:
        config["results"] + "gvcf/{batch}/{seq}{indiv}_{library}.chr{chr}.log",
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
        vcfs = lambda wildcards: expand(config["results"] + "gvcf/{batch}/{seq}{indiv}_{library}.chr{chr}.g.vcf.gz",
            chr=CHROMOSOMES,
            batch=wildcards.batch,
            seq = wildcards.seq,
            indiv = wildcards.indiv,
            library = wildcards.library),
    output:
        vcf = config["results"] + "gvcf/{batch}/{seq}{indiv}_{library}.g.vcf.gz",
        #tbi = config["results"] + "gvcf/{seq}{indiv}_{library}.g.vcf.gz.tbi",
    conda: "../envs/common.yaml"
    threads: 1
    resources: nodes = 1
    shell: """
        bcftools concat {input.vcfs} \
            -o {output.vcf} \
            -Oz \
        """

# Consolidation of GVCFs
rule create_sample_map:
    """Create sample map that contains names and paths to all VCFs to be used in consolidate rule."""
    input:
        gvcfs = lambda wildcards: expand("{results}gvcf/{sample}.chr{chr}.g.vcf.gz",
            results=config["results"],
            sample=collect_samples(fmt="{batch}/{seq}{indiv}_{library}"),
            chr=wildcards.chr),
    output:
        sample_map = temp(config["results"] + "db/sample-maps/{dataset}.chr{chr}.sample-map"),
    params:
        json = lambda wildcards: config["results"] + f"db/{wildcards.dataset}/{wildcards.chr}/callset.json",
    threads: 1
    resources: nodes = 1
    run:
        import json
        print(f'{params.json}')
        try:
            with open(f'{params.json}') as f:
                print('try beginning')
                data = json.load(f)
                samples_already_in_datastore = [sample['sample_name'] for sample in data['callsets']]
                print('try end')
        except FileNotFoundError:
            print('except')
            samples_already_in_datastore = []
        print(samples_already_in_datastore)

        with open(output.sample_map, "w") as sample_map:
            for gvcf in input.gvcfs:
                sample = gvcf.split("/")[-1].split(".")[0]
                print(sample)
                if sample not in samples_already_in_datastore:
                    print("not already in datastore")
                    sample_map.write(f"{sample}\t{gvcf}\n")

rule consolidate:
    """Combine the chromosomes of .g.vcf files into GenomicsDB datastore."""
    input:
        gvcfs = lambda wildcards: expand("{results}gvcf/{sample}.chr{chr}.g.vcf.gz{ext}",
        #gvcfs = lambda wildcards: expand("{results}gvcf/{sample}.g.vcf.gz{ext}",
            results=config["results"],
            sample=collect_samples(fmt="{batch}/{seq}{indiv}_{library}"),
            chr=wildcards.chr,
            ext=["", ".tbi"],),
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
rule joint_call_cohort:
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
        split = config["results"] + "joint_call/biallelic/{dataset}.{mode}.chr{chr}.bcf",
    params:
        #-e'type{params.equality}"snp"' \
        #equality = lambda wildcards: "=" if wildcards.mode == "indel" else "!=",
        equality_option = lambda wildcards: """-e'type="snp"'""" if wildcards.mode == "indel" else (
            """-e'type!="snp"'""" if wildcards.mode == "SNP" else ""
            ),
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
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
            -Ob \
            -o {output.split} \
        """
