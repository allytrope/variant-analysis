"""Contain rules for genotype refinement and haplotype estimation.
This comes after hard_filter.smk or can be modified to come after variant_recalibration.smk."""

# -------------
# Use either of the following two rules

rule deduplicate_individuals:
    """Keep one of each animal and then simply names down to ids.
    Samples are prioritizes in the order: WGS, WES, GBS, AMP.
    Then the most by the last run name alphabetically (which is usually the newest)."""
    input:
        #bcf = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.SNP.chr{chr}.bcf",
        vcf = config["results"] + "hard_filtered/pass/{dataset}.{mode}.chr{chr}.bcf",
        tbi = config["results"] + "hard_filtered/pass/{dataset}.{mode}.chr{chr}.bcf.csi",
        # vcf = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.vcf.gz",
        # tbi = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
    output:
        vcf = config["results"] + "genotypes/deduplicated/{dataset}.{mode}.chr{chr}.vcf.gz",
        tmp_samples = temp(config["results"] + "genotypes/deduplicated/{dataset}.{mode}.chr{chr}.samples.list"),
        tmp_samples2 = temp(config["results"] + "genotypes/deduplicated/{dataset}.{mode}.chr{chr}.samples2.list"),
    conda: "../envs/common.yaml"
    # Split full sample names into seq + animal_id + run_id assuming a field like so: {seq}{animal_id}_{run_id}
    # Sort by priority
    # Only keep first of each animal_id
    # Then apply list to VCF
    # Reheader sample names to individual names
    shell: """
        bcftools query -l {input.vcf} \
            | sed 's/.\{{3\}}/&\t/; s/_/\t/' \
            | csvtk sort \
                -k 2,1:u,3:r \
                -L 1:<(for i in WGS lpWGS WES GBS AMP; do echo $i; done) \
                -t \
                -H \
            | awk 'BEGIN {{FS="\t"; OFS=""}} $2!=prev {{print $1,$2,"_"$3}} {{prev=$2}}' \
            | tee {output.tmp_samples} \
            | cut -d '_' -f 1 \
            | cut -c 4- \
            > {output.tmp_samples2}; \
        bcftools view {input.vcf} \
            -S {output.tmp_samples} \
        | bcftools reheader  \
            -s {output.tmp_samples2} \
        | bcftools view \
            -Oz \
        > {output.vcf} \
        """

rule only_one_seq_type:
    """Subset down to only one sample for each individual and of only a specific sequencing type."""
    input:
        vcf = config["results"] + "hard_filtered/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "hard_filtered/pass/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
    output:
        vcf = config["results"] + "genotypes/only_{seq}/{dataset}.{mode}.chr{chr}.vcf.gz",
        tmp_samples = temp(config["results"] + "genotypes/only_{seq}/{dataset}.{mode}.chr{chr}.samples.list"),
        tmp_samples2 = temp(config["results"] + "genotypes/only_{seq}/{dataset}.{mode}.chr{chr}.samples2.list"),
    conda: "../envs/common.yaml"
    shell: """
        bcftools query -l {input.vcf} \
            | sed 's/.\{{3\}}/&\t/; s/_/\t/' \
            | csvtk grep -p {wildcards.seq} -t -f 1 \
            | awk 'BEGIN {{FS="\t"; OFS=""}} $2!=prev {{print $1,$2,"_"$3}} {{prev=$2}}' \
            | tee {output.tmp_samples} \
            | cut -d '_' -f 1 \
            | cut -c 4- \
            > {output.tmp_samples2}; \
        bcftools view {input.vcf} \
            -S {output.tmp_samples} \
        | bcftools reheader  \
            -s {output.tmp_samples2} \
        | bcftools view \
            --min-ac 1 \
            -Oz \
        > {output.vcf} \
        """

# ----------------

rule genotype_posteriors:
    """Calculate genotype posterior probabilties. This adds a PP field in the FORMAT column of the VCF for the genotype posteriors.
    Also, the GQ and GT fields may also be subject to change. These are calculated using information from the other samples in the file
    and also optionally from trio information."""
    # Switch input vcf/tbi as needed
    input:
        # TODO: Generalize the "only_{seq}"
        #vcf = config["results"] + "genotypes/only_WES/{dataset}.{mode}.chr{chr}.vcf.gz",
        #tbi = config["results"] + "genotypes/only_WES/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
        vcf = config["results"] + "genotypes/deduplicated/{dataset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "genotypes/deduplicated/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
        #ped = config["resources"] + "pedigree/trios.ped",
        #ped = config["resources"] + "pedigree/plink.ped",
    output:
        config["results"] + "genotypes/posteriors/{dataset}.{mode}.chr{chr}.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    # --tmp-dir ~/tmp/{rule}
    #--pedigree {input.ped} \
    shell: """
        gatk --java-options "-Xmx8g" CalculateGenotypePosteriors \
            -V {input.vcf} \
            -O {output} \
            """

rule genotype_filtration:
    """Filter genotypes by GQ.
    This adds the filter tag if fails and sets "./." as new genotype."""
    input:
        vcf = config["results"] + "genotypes/posteriors/{dataset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "genotypes/posteriors/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
    output:
        vcf = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.vcf.gz",
    params:
        GQ = config["filtering"]["GQ"],
        DP = 5,
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    shell: """
        gatk --java-options "-Xmx8g" VariantFiltration \
            -V {input.vcf} \
            --genotype-filter-name "GQ{params.GQ}" \
            --genotype-filter-expression "GQ < {params.GQ}" \
            --genotype-filter-name "DP{params.DP}" \
            --genotype-filter-expression "DP < {params.DP}" \
            --set-filtered-genotype-to-no-call \
            -O {output.vcf}"""

rule genotype_passing:
    """Remove variants that don't have a low alternate allele frequency.
    A min_AC of 1 is necessary to remove ACs that were set to 0 during genotype refinement."""
    # wildcard_constraints:
    #     #subset = "[^founders2_MAF1|common_between_founding_cohorts2|founding_cohorts1_1|founding_cohorts1_2]"
    #     subset = "all2"
    input:
        #vcf = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.vcf.gz",
        vcf = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.vcf.gz",
        #subpop = config["subpop"],
        subpop = lambda wildcards: branch(
            True if "all" not in wildcards.subset else False,
            #then = config["subpop"],
            then = config["resources"] + f"subpop/{wildcards.subset}.list",
        )
    output:
        bcf = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.bcf",
    params:
        #samples = ','.join(SAMPLES),
        min_AF = 0, #config["filtering"]["min_AF"],
        max_AF = 1, #config["filtering"]["max_AF"],
        min_AC = 1, #config["filtering"]["min_AC"],  # Must be 1 or greater
        subset_samples = lambda wildcards, input: "" if "all" in wildcards.subset else f"-S {input.subpop}",
    threads: 1
    resources: nodes = 1

    # -e 'F_MISSING > 0.1' \
    # shell: """
    #     bcftools view {input.vcf} \
    #         {params.subset_samples} \
    #         -Ou \
    #     | bcftools view \
    #         --min-af {params.min_AF} \
    #         --max-af {params.max_AF} \
    #         --min-ac {params.min_AC} \
    #         -e 'F_MISSING>0.2' \
    #         -Ob \
    #         -o {output.bcf} \
    #     """
    shell: """
        bcftools view {input.vcf} \
            {params.subset_samples} \
            -Ou \
        | bcftools view \
            --min-af {params.min_AF} \
            --max-af {params.max_AF} \
            --min-ac {params.min_AC} \
            -Ob \
            -o {output.bcf} \
        """

rule prune_LD:
    """Remove SNVs with high LD."""
    input:
        bcf = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.bcf",
    output:
        bcf = config["results"] + "genotypes/pruned/{dataset}.{subset}.{mode}.chr{chr}.bcf",
    threads: 1
    resources:
        nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools +prune {input.bcf} \
            -m 0.6 \
            --random-seed 7340 \
            -Ob \
            -o {output.bcf} \
        """

rule create_plink_files:
   """Create PLINK .bed, .bim, and .fam files."""
    input:
        bcf = "{path}/{dataset}.{subset}.{mode}.chr{chr}.bcf",
        demographics = config["resources"] + "pedigree/demographics.tsv",
    output:
        bed = "{path}/plink/{dataset}.{subset}.{mode}.chr{chr}.bed",
        bim = "{path}/plink/{dataset}.{subset}.{mode}.chr{chr}.bim",
        fam = "{path}/plink/{dataset}.{subset}.{mode}.chr{chr}.fam",
        update_parents = "{path}/plink/{dataset}.{subset}.{mode}.chr{chr}.update_parents.tsv",
        update_sex = "{path}/plink/{dataset}.{subset}.{mode}.chr{chr}.update_sex.tsv",
    params:
        out_prefix = subpath(output.bed, strip_suffix=".bed"),
    conda: "../envs/rvtests.yaml"
    # --chr {wildcards.chr} \
    shell: """
        cat {input.demographics} \
        | csvtk cut -f Id,Sire,Dam -t \
        | awk 'BEGIN {{FS="\t"; OFS="\t"}} NR!=1 {{print "F1",$0}} NR==1 {{print "Family",$0}}' \
        | sed 's/\t\t/\t0\t/g' \
        | sed 's/\t$/\t0/g' \
        > {output.update_parents}; \

        cat {input.demographics} \
        | csvtk cut -f Id,Sex -t \
        | awk 'BEGIN {{FS="\t"; OFS="\t"}} NR!=1 {{print "F1",$0}} NR==1 {{print "Family",$0}}' \
        | sed 's/Male/1/g;s/Female/2/g;s/Unknown/0/g' \
        > {output.update_sex}; \

        plink \
            --bcf {input.bcf} \
            --const-fid F1 \
            --recode \
            --make-bed \
            --out {params.out_prefix} \
            --update-parents {output.update_parents} \
            --update-sex {output.update_sex} \
            --allow-extra-chr \
        """
        


# Minor differences
rule create_autosomal_plink_files:
   """Create PLINK .bed, .bim, and .fam files."""
    input:
        bcf = "{path}/{dataset}.{subset}.{mode}.autosomal.bcf",
        demographics = config["resources"] + "pedigree/demographics.tsv",
    output:
        bed = "{path}/plink/{dataset}.{subset}.{mode}.autosomal.bed",
        bim = "{path}/plink/{dataset}.{subset}.{mode}.autosomal.bim",
        fam = "{path}/plink/{dataset}.{subset}.{mode}.autosomal.fam",
        update_parents = "{path}/plink/{dataset}.{subset}.{mode}.autosomal.update_parents.tsv",
        update_sex = "{path}/plink/{dataset}.{subset}.{mode}.autosomal.update_sex.tsv",
    params:
        out_prefix = subpath(output.bed, strip_suffix=".bed"),
    conda: "../envs/rvtests.yaml"
    shell: """
        cat {input.demographics} \
        | csvtk cut -f Id,Sire,Dam -t \
        | awk 'BEGIN {{FS="\t"; OFS="\t"}} NR!=1 {{print "F1",$0}} NR==1 {{print "Family",$0}}' \
        | sed 's/\t\t/\t0\t/g' \
        | sed 's/\t$/\t0/g' \
        > {output.update_parents}; \

        cat {input.demographics} \
        | csvtk cut -f Id,Sex -t \
        | awk 'BEGIN {{FS="\t"; OFS="\t"}} NR!=1 {{print "F1",$0}} NR==1 {{print "Family",$0}}' \
        | sed 's/Male/1/g;s/Female/2/g;s/Unknown/0/g' \
        > {output.update_sex}; \

        plink \
            --bcf {input.bcf} \
            --const-fid F1 \
            --recode \
            --make-bed \
            --out {params.out_prefix} \
            --update-parents {output.update_parents} \
            --update-sex {output.update_sex} \
            --allow-extra-chr \
        """

















# rule exonic_regions_only_of_genotype_passing:
#     """Subset to WES regions determined by mosdepth from genotype passing VCF."""
#     input:
#         bcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.bcf",
#         bed = config["results"] + "coverage/common_WES_0.5_loci.bed",  # TODO: Generalize this
#     output:
#         bcf = config["results"] + "genotypes/pass/exonic_regions/{dataset}.{mode}.chr{chr}.bcf",
#     shell: """
#         bcftools view {input.bcf} \
#             -R {input.bed} \
#             -Oz \
#             -o {output.bcf} \
#         """

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
    conda: "../envs/common.yaml"
    shell: """
        bcftools annotate -x INFO,^FORMAT/GT {input.trio} \
        | bcftools view -H -i "F_MISSING==0" \
        | awk '{{print $10,$11,$12}}' \
        | python3 {params.script} \
        | sort \
        | uniq -c > {output.counts} \
        """

rule whatshap_cohort:
    """Phase GATK VCF. Uses BAMs and parent information to phase along parent genomes."""
    input:
        vcf = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.bcf",
        idx = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.bcf.csi",
        ref_fasta = config["ref_fasta"],
        ref_idx = config["ref_fasta"] + ".fai",
        # Child, sire, and dam BAMs and BAM idxs also required, but specified under "params" and "shell". These should already exist from earlier commands anyway.
    output:
        phased = temp(config["results"] + "haplotypes/whatshap/trios/chr{chr}/{dataset}.{sample}.{mode}.chr{chr}.vcf.gz"),
    params:
        child_bam = lambda wildcards, input: find_bam(0, wildcards, input),
        sire_bam = lambda wildcards, input: find_bam(1, wildcards, input),
        dam_bam = lambda wildcards, input: find_bam(2, wildcards, input),
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        whatshap phase {input.vcf} {params.child_bam} {params.sire_bam} {params.dam_bam} \
            --reference={input.ref_fasta} \
            --tag=PS \
            -o {output.phased} \
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
        # vcf = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz",
        # idx = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz.tbi",
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
    conda: "../envs/common.yaml"
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
    conda: "../envs/common.yaml"
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
    conda: "../envs/common.yaml"
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
    conda: "../envs/common.yaml"
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