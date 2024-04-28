"""Contain rules for genotype refinement and haplotype estimation.
This comes after hard_filter.smk or can be modified to come after variant_recalibration.smk."""


rule genotype_posteriors:
    """Calculate genotype posterior probabilties. This adds a PP field in the FORMAT column of the VCF for the genotype posteriors.
    Also, the GQ and GT fields may also be subject to change. These are caluclated using information from the other samples in the file
    and also optionally from trio information."""
    input:
        vcf = config["results"] + "hard_filtered/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "hard_filtered/pass/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
        #ped = config["resources"] + "pedigree/trios.ped",
    output:
        config["results"] + "genotypes/posteriors/{dataset}.{mode}.chr{chr}.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    # --tmp-dir ~/tmp/{rule}
    shell: """
        gatk --java-options "-Xmx8g" CalculateGenotypePosteriors \
            -V {input.vcf} \
            -O {output} \
            """
            #-ped {input.ped} \

rule genotype_filtration:
    """Filter genotypes by GQ.
    This adds the filter tag if fails and sets "./." as new genotype."""
    input:
        vcf = config["results"] + "genotypes/posteriors/{dataset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "genotypes/posteriors/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
    output:
        vcf = temp(config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.vcf.gz"),
    params:
        GQ = config["filtering"]["GQ"],
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    shell: """
        gatk --java-options "-Xmx8g" VariantFiltration \
            -V {input.vcf} \
            --genotype-filter-expression "GQ < {params.GQ}" \
            --genotype-filter-name "GQ{params.GQ}" \
            --set-filtered-genotype-to-no-call \
            -O {output.vcf}"""

rule genotype_passing:
    """Remove variants that don't have a low alternate alleles.
    A min_AC of 1 is necessary to remove ACs that were set to 0 during genotype refinement."""
    input:
        vcf = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.vcf.gz",
    output:
        bcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.bcf",
    params:
        samples = ','.join(SAMPLES),
        min_AF = config["filtering"]["min_AF"],
        max_AF = config["filtering"]["max_AF"],
        min_AC = config["filtering"]["min_AC"],  # Must be 1 or greater
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools view {input.vcf} \
            -s {params.samples} \
            -Ou \
        | bcftools view \
            --min-af {params.min_AF} \
            --max-af {params.max_AF} \
            --min-ac {params.min_AC} \
            -Ob \
            -o {output.bcf} \
        """

rule count_rates_of_Mendelian_errors_by_GQ:
    """Find rate of Mendelian errors. For quality checking only."""
    input:
        trio = config["results"] + "genotypes/pass/trios/{dataset}.{sample}.{mode}.chr{chr}.vcf.gz",
    output:
        counts = config["results"] + "genotypes/pass/counts/{dataset}.{sample}.{mode}.chr{chr}.counts",
    params:
        script = "workflow/scripts/count_Mendelian_errors.py",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools annotate -x INFO,^FORMAT/GT {input.trio} \
        | bcftools view -H -i "F_MISSING==0" \
        | awk '{{print $10,$11,$12}}' \
        | python3 {params.script} \
        | sort \
        | uniq -c > {output.counts} \
        """

def find_bam(idx, wildcards, input):
    """Find BAM file from sample name."""
    with open(input.tsv, "r") as f:
        for line in f.readlines():
            if line.startswith(wildcards.sample):
                trio = line.strip("\n").split("\t")
    sample = trio[idx]
    if sample == "":
        return ""
    else:
        sample_runs = []
        for run in SAMPLE_RUNS:
            if sample in run:
                sample_runs.append(config["results"] + "alignments/recalibrated/" + run + ".bam")
        return sample_runs
rule whatshap_trio:
    """Haplotype assembly to be used as phase set for SHAPEIT4. Uses BAMs and parent information to phase along parent genomes."""
    input:
        #vcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
        #idx = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
        vcf = config["results"] + "genotypes/pass/trios/chr{chr}/{dataset}.{sample}.{mode}.chr{chr}.vcf.gz",
        idx = config["results"] + "genotypes/pass/trios/chr{chr}/{dataset}.{sample}.{mode}.chr{chr}.vcf.gz.tbi",
        ref_fasta = config["ref_fasta"],
        ref_idx = config["ref_fasta"] + ".fai",
        tsv = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.tsv",
        ped = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.ped",
        # Child, sire, and dam BAMs and BAM idxs also required, but specified under "params" and "shell". These should already exist from earlier commands anyway.
    output:
        phased = temp(config["results"] + "haplotypes/whatshap/trios/chr{chr}/{dataset}.{sample}.{mode}.chr{chr}.vcf.gz"),
    params:
        child_bam = lambda wildcards, input: find_bam(0, wildcards, input),
        sire_bam = lambda wildcards, input: find_bam(1, wildcards, input),
        dam_bam = lambda wildcards, input: find_bam(2, wildcards, input),
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        whatshap phase {input.vcf} {params.child_bam} {params.sire_bam} {params.dam_bam} \
            --reference={input.ref_fasta} \
            --ped={input.ped} \
            --tag=PS \
            -o {output.phased} \
        """

rule indiv_vcf:
    """Take only child from trio VCF."""
    input:
        trio = config["results"] + "haplotypes/whatshap/trios/chr{chr}/{dataset}.{sample}.{mode}.chr{chr}.vcf.gz",
    output:
        indiv = config["results"] + "haplotypes/whatshap/indivs/chr{chr}/{dataset}.{sample}.{mode}.chr{chr}.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        SAMPLE=$(bcftools query -l {input.trio} | head -n 1); \
        bcftools view {input.trio} \
            -s $SAMPLE \
            -Ou \
            -o {output.indiv} \
        """

rule merge_whatshap:
    input:
        trios = lambda wildcards: expand(config["results"] + "haplotypes/whatshap/indivs/chr{chr}/{dataset}.{sample}.{mode}.chr{chr}.vcf.gz",
            dataset=wildcards.dataset,
            sample=SAMPLES,
            mode=wildcards.mode,
            chr=wildcards.chr),
        idxs = lambda wildcards: expand(config["results"] + "haplotypes/whatshap/indivs/chr{chr}/{dataset}.{sample}.{mode}.chr{chr}.vcf.gz.tbi",
            dataset=wildcards.dataset,
            sample=SAMPLES,
            mode=wildcards.mode,
            chr=wildcards.chr),
    output:
        merged = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools merge {input.trios} \
            -Oz \
            -o {output.merged} \
        """

rule restrict_vcf_seq_type:
    """Create VCF with only one type of sequencing method (e.g., WES, WGS)"""
    input:
        vcf = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz",
    output:
        samples = temp(config["results"] + "haplotypes/whatshap/{dataset}.{mode}.chr{chr}.{seq}_samples.list"),
        vcf = config["results"] + "haplotypes/whatshap/{seq}/{dataset}.{mode}.chr{chr}.vcf.gz",
    conda: "../envs/bio.yaml"
    shell: """
        bcftools query -l {input.vcf} | grep {wildcards.seq} > {output.samples};
        bcftools view {input.vcf} \
            -S {output.samples} \
            -Oz \
            -o {output.vcf} \
        """

rule WES_only:
    """Subset out WES samples and keep WES regions determined by mosdepth."""
    input:
        vcf = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz",
        bed = config["results"] + "coverage/mosdepth/common_WES_0.5_loci.bed",  # TODO: Generalize this
    output:
        samples = temp(config["results"] + "haplotypes/whatshap/{dataset}.{mode}.chr{chr}.{seq}_samples.list"),
        vcf = config["results"] + "haplotypes/whatshap/{seq}_SNPs/{dataset}.{mode}.chr{chr}.vcf.gz",
    shell: """
        bcftools query -l {input.vcf} | grep {wildcards.seq} > {output.samples};
        bcftools view {input.vcf} \
            -S {output.samples} \
            -R {input.bed} \
            -Oz \
            -o {output.vcf} \
        """