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
    shell: """
        gatk --java-options "-Xmx8g" CalculateGenotypePosteriors \
            -V {input.vcf} \
            -O {output} \
            --tmp-dir ~/tmp/{rule}"""
            #-ped {input.ped} \

rule genotype_filtration:
    """Filter genotypes by GQ.
    This adds the filter tag if fails and sets "./." as new genotype."""
    input:
        vcf = config["results"] + "genotypes/posteriors/{dataset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "genotypes/posteriors/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
    output:
        vcf = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.vcf.gz",
    params:
        GQ = 45,
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
        samples = config["samples"],
    output:
        vcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
    params:
        min_AF = config["filtering"]["min_AF"],
        max_AF = config["filtering"]["max_AF"],
        min_AC = config["filtering"]["min_AC"],  # Should be 1 or greater
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools view {input.vcf} \
            -S {input.samples} \
            -Ou \
        | bcftools view \
            --min-af {params.min_AF} \
            --max-af {params.max_AF} \
            --min-ac {params.min_AC} \
            -Oz \
            -o {output.vcf} \
        """

rule largest_seq_per_organism:
    """Create a list of samples where only one sample from each organism is kept.
    Only the sequencing type with the most data is kept, which is determined by WGS > WES > GBS > AMP.
    This ordering happens to work here because of the alphabetical ordering.
    Used for determining which sample to use as parent in later step."""
    input:
        samples = config["samples"],
    output:
        largest_samples = config["results"] + "haplotypes/pedigree/{dataset}.largest_samples.list",
    # sed separates id and sequence type. E.g. WGS12345 -> WGS    12345
    # awk flips columns
    # sort rows
    # merge rows based on first column
    # take largest id with largest sequence type (this works because "GBS", "WES", "WGS" are in alphabetical order)
    # grep to remove sample names with an underscore
    # sort
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        sed 's/./&\t/3' {input.samples} \
        | awk -v OFS='\t' '{{print $2,$1}}' \
        | sort \
        | awk '$1!=p{{if(p)print s; p=$1; s=$0; next}}{{sub(p,x); s=s $0}} END {{print s}}' \
        | awk -v OFS='' '{{print $NF,$1}}' \
        | grep -v "_" \
        | sort > {output.largest_samples} \
        """

rule add_seq_to_children:
    """Add seq type (AMP, GBS, WES, and/or WGS) to individual ids in pedigree
    and repeat entries if there are multiple sequencing types for an individual."""
    input:
        samples = config["samples"],
        parents = config["pedigree"],
        largest_samples = config["results"] + "haplotypes/pedigree/{dataset}.largest_samples.list",
    output:
        parents = temp(config["results"] + "haplotypes/pedigree/{dataset}.children_with_seq.tsv"),
    threads: 1
    resources: nodes = 1
    shell: """
        sed 's/AMP/AMP\t/g;s/GBS/GBS\t/g; s/WES/WES\t/g; s/WGS/WGS\t/g' {input.samples} \
        | sort -k 2 \
        | join - {input.parents} -1 2 -2 1 \
        | awk '{{print $2$1"\t"$3"\t"$4}}' \
        > {output.parents} \
        """
    
rule add_seq_to_parents:
    """Prepend sequence type to parents in pedigree.
    Uses WGS if available, otherise tries the same sequence type as child.
    If the parent is not sequenced, no sequence type is prepended."""
    input:
        samples = config["samples"],
        parents = config["results"] + "haplotypes/pedigree/{dataset}.children_with_seq.tsv",
    output:
        tsv = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.tsv",
    threads: 1
    resources: nodes = 1
    run:
        samples = None
        with open(input.samples, "r") as f:
            samples = f.read().strip().split("\n")
        with open(input.parents, "r") as f:
            with open(output.tsv, "a") as out:
                for line in f:
                    child, sire, dam = line.strip("\n").split("\t")
                    seq_type = child[:3]
                    if f"WGS{sire}" in samples:
                        sire = f"WGS{sire}"
                    elif f"{seq_type}{sire}" in samples:
                        sire = f"{seq_type}{sire}"
                    else:
                        sire = ""
                    if f"WGS{dam}" in samples:
                        dam = f"WGS{dam}"
                    elif f"{seq_type}{dam}" in samples:
                        dam = f"{seq_type}{dam}"
                    else:
                        dam = ""
                    out.write(f"{child}\t{sire}\t{dam}\n")

rule make_forced_ped_format:
    """Add fields to make the correct number of columns for a PLINK PED file.
    However, these extra fields don't actually hold any relevant information.
    They are just a requirement for WhatsHap."""
    input:
        tsv = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.tsv",
    output:
        ped = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.ped",
    threads: 1
    resources: nodes = 1
    shell: """
        awk 'BEGIN {{OFS="\t"}} {{print 0,$1,$2,$3,0,0}}' {input.tsv} \
        | sed 's/\t\t/\t0\t/g' \
        | sed 's/\t\t/\t0\t/g' \
        > {output.ped} \
        """

rule trios_only_tsv:
    input:
        tsv = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.tsv",
    output:
        tsv = config["results"] + "haplotypes/pedigree/{dataset}.trios_only.tsv",
    threads: 1
    resources: nodes = 1
    shell: """
        grep -E "(WGS.*|WES.*|GBS.*|AMP.*){3}" {input.tsv} \
        > {output.tsv} \
        """

rule make_trio_vcf:
    """Keep only one child-sire-dam trio."""
    input:
        vcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
        tsv = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.tsv",
    output:
        trio = config["results"] + "genotypes/pass/trios/chr{chr}/{dataset}.{sample}.{mode}.chr{chr}.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        read -r CHILD SIRE DAM <<< $(grep ^{wildcards.sample} {input.tsv}); \
        if [ -n "$SIRE" ]; \
        then SIRE=",$SIRE"; \
        fi; \
        if [ -n "$DAM" ]; \
        then DAM=",$DAM"; \
        fi; \
        bcftools view {input.vcf} \
            -s {wildcards.sample}$SIRE$DAM \
            -Oz \
            -o {output.trio} \
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
        return config["results"] + f"alignments/recalibrated/{sample}.bam"
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
        phased = config["results"] + "haplotypes/whatshap/trios/chr{chr}/{dataset}.{sample}.{mode}.chr{chr}.vcf.gz",
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
            sample=SAMPLE_NAMES,
            mode=wildcards.mode,
            chr=wildcards.chr),
        idxs = lambda wildcards: expand(config["results"] + "haplotypes/whatshap/indivs/chr{chr}/{dataset}.{sample}.{mode}.chr{chr}.vcf.gz.tbi",
            dataset=wildcards.dataset,
            sample=SAMPLE_NAMES,
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

rule make_scaffold:
    """Create scaffold for SHAPEIT4."""
    input:
        vcf = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
        fam = config["results"] + "haplotypes/pedigree/all_samples.trios_only.tsv",
    output:
        scaffold = config["results"] + "haplotypes/scaffolds/{dataset}.{mode}.chr{chr}.vcf.gz",
    threads: 1
    resources: nodes = 1
    shell: """
        makeScaffold \
            --gen {input.vcf} \
            --fam {input.fam} \
            --reg {wildcards.chr} \
            --out {output.scaffold} \
        """

rule shapeit4_imputation:
    """Haplotype estimation and imputation."""
    input:
        vcf = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz",
        csi = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz.csi",
        scaffold = config["results"] + "haplotypes/scaffolds/{dataset}.{mode}.chr{chr}.vcf.gz",
        scaffold_csi = config["results"] + "haplotypes/scaffolds/{dataset}.{mode}.chr{chr}.vcf.gz.csi",
    output:
        phased = config["results"] + "haplotypes/SHAPEIT4/with_scaffold/{dataset}.{mode}.chr{chr}.vcf.gz",
    # When using SHAPIT4.1 or greater, --map not required (though surely still helpful).
    # chr appears to still be mandatory even if there is only one chromosome in file.
    params:
        PS = 0.0001,  # Recommended value by SHAPEIT4
    log: config["results"] + "haplotypes/SHAPEIT4/log/{dataset}.{mode}.chr{chr}.log",
    threads: 1
    resources: nodes = 1
    conda: "../envs/shapeit4.yaml"
    shell: """
        shapeit4 \
            --input {input.vcf} \
            --scaffold {input.scaffold} \
            --region {wildcards.chr} \
            --sequencing \
            --use-PS {params.PS} \
            --output {output.phased} \
            --log {log} \
            --thread {threads} \
        """

rule add_annotations:
    """Adds FORMAT annotations that were removed during processing with SHAPEIT4."""
    input: 
        vcf = config["results"] + "haplotypes/whatshap/all/{dataset}.{mode}.chr{chr}.vcf.gz",
        phased = config["results"] + "haplotypes/SHAPEIT4/{dataset}.{mode}.chr{chr}.vcf.gz",
    output:
        annotated = config["results"] + "haplotypes/SHAPEIT4/annotated/{dataset}.{mode}.chr{chr}.vcf.gz",
    threads: 1 
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools annotate {input.vcf} \
            -a {input.phased} \
            -c FORMAT/GT \
            -o {output.annotated} \
            -Oz \
        """
