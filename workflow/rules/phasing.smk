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
    """Filter genotypes for GQ < 20 (keeping only genotypes with >99% chance of being correct).
    This adds the filter tag if fails and sets "./." as new genotype."""
    input:
        vcf = config["results"] + "genotypes/posteriors/{dataset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "genotypes/posteriors/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
    output:
        vcf = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/gatk.yaml"
    shell: """
        gatk --java-options "-Xmx8g" VariantFiltration \
            -V {input.vcf} \
            --genotype-filter-expression "GQ < 20" \
            --genotype-filter-name "GQ20" \
            --set-filtered-genotype-to-no-call \
            -O {output.vcf}"""

rule only_SNPs:
    """Get rid of rows with AC=0."""
    input:
        vcf = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.vcf.gz",
    output:
        vcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
    shell: """
        bcftools view {input.vcf} \
            --min-ac 1 \
            -Oz \
            -o {output.vcf}
        """

rule AC_filter_and_subset_samples:
    """Remove variants that don't have a low alternate alleles."""
    input:
        vcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
        samples = config["samples"],
    output:
        vcf = config["results"] + "genotypes/subsets/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz",
    params:
        min_AF = 0.05,
        max_AF = 0.95,
        min_AC = 10,
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

# This rule works when run manually, but not when run as a rule.
rule create_subset_list:
    """Create a list of samples where only one sample from each organism is kept.
    Only the sequencing type with the most data is kept, which is determined by WGS > WES > GBS."""
    input:
        #all_samples = config["samples"],
        vcf = config["results"] + "genotypes/filtered/{dataset}.SNP.chr1.filtered.vcf.gz",
    output:
        largest_samples = config["results"] + "genotypes/subsets/{dataset}.largest_samples.list",
    # bcftools retrieves samples
    # sed separates id and sequence type. E.g. WGS12345 -> WGS    12345
    # awk flips columns
    # sort rows
    # merge rows based on first column
    # take largest id with largest sequence type (this works because "GBS", "WES", "WGS" are in alphabetical order)
    # grep to remove sample names with an underscore
    # sort
    conda: "../envs/bio.yaml"
    shell: """
        bcftools query -l {input.vcf} \
            | sed 's/./&\t/3' \
            | awk -v OFS='\t' '{{print $2,$1}}' \
            | sort \
            | awk '$1!=p{{if(p)print s; p=$1; s=$0; next}}{{sub(p,x); s=s $0}} END {{print s}}' \
            | awk -v OFS='' '{{print $NF,$1}}' \
            | grep -v "_" \
            | sort > {output.largest_samples} \
        """

### Phasing

# Not currently being used
rule whatshap:
    """Haplotype assembly to be used as phase set for SHAPEIT4."""
    input: 
        vcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
        bams = expand(config["results"] + "alignments_recalibrated/{sample}.bam", sample=SAMPLE_NAMES),  # Should include all samples in VCF
        ref_fasta = config["ref_fasta"],
        ped = ""  # PLINK PED file (but only uses columns 2, 3, and 4)
    output:
        phased = config["results"] + "haplotypes/whatshap/{dataset}.{mode}.chr{chr}.phased.vcf.gz",
    threads: 1
    resources: nodes = 1
    shell: """
        whatshap phase {input.vcf} {input.bams} \
            --reference={input.ref_fasta} \
            --ped {input.ped} \
            -o {output.vcf} \
        """

rule add_labels_to_pedigree:
    """Replace animal ids with those including sequence type. For examples, 12345 -> GBS12345 or 11111 -> WGS11111."""
    input:
        vcf = config["results"] + "genotypes/subsets/GBS_WES_WGS_labels.284_and_5_swapped.SNP.chr11.min_AF0.01.min_AC2.vcf.gz",
        pedigree = config["resources"] + "pedigree/{pedigree}.txt",
    output:
        pedigree = config["resources"] + "pedigree/{pedigree}.labels.txt",
    run:
        samples = shell(f"bcftools query -l {input.vcf}", iterable=True)
        with open(input.pedigree) as f:
            text = f.read()
        for sample in samples:
            text = text.replace(sample[3:], sample)
        with open(output.pedigree, "x") as f:
            f.write(text)

rule make_scaffold:
    """Create scaffold for SHAPEIT4."""
    input:
        #vcf = config["results"] + "genotypes/filtered/{dataset}.SNP.chr{chr}.filtered.vcf.gz",  # Taking hard filtered .vcf
        #tbi = config["results"] + "genotypes/filtered/{dataset}.SNP.chr{chr}.filtered.vcf.gz.tbi",
        vcf = config["results"] + "genotypes/subsets/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz",
        tbi = config["results"] + "genotypes/subsets/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz.tbi",
        fam = config["resources"] + "pedigree/parents.labels.txt",
    output:
        scaffold = config["results"] + "haplotypes/scaffolds/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz",
    threads: 1
    resources: nodes = 1
    shell: """
        makeScaffold \
            --gen {input.vcf} \
            --fam {input.fam} \
            --reg {wildcards.chr} \
            --out {output.scaffold} \
        """

rule shapeit4_without_ref:
    """Haplotype estimation and imputation. Not using reference this time, but still using scaffold."""
    input:
        #vcf = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.filtered.vcf.gz",  # Taking hard filtered .vcf
        #csi = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.filtered.vcf.gz.csi",
        vcf = config["results"] + "genotypes/subsets/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz",
        csi = config["results"] + "genotypes/subsets/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz.csi",
        scaffold = config["results"] + "haplotypes/scaffolds/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz",
        scaffold_csi = config["results"] + "haplotypes/scaffolds/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz.csi",
    output:
        phased = config["results"] + "haplotypes/SHAPEIT4/{dataset}.{subset}.{mode}.chr{chr}.no_ref.vcf.gz",
    # When using SHAPIT4.1 or greater, --map not required (though surely still helpful).
    # chr appears to still be mandatory even if there is only one chromosome in file.
    # --scaffold is useful for incorporation pedigree information.
    log: config["results"] + "haplotypes/SHAPEIT4/log/{dataset}.{subset}.{mode}.chr{chr}.no_ref.log",
    threads: 24
    resources: nodes = 24
    conda: "../envs/shapeit4.yaml"
    shell: """
        shapeit4 \
            --input {input.vcf} \
            --scaffold {input.scaffold} \
            --region {wildcards.chr} \
            --sequencing \
            --output {output.phased} \
            --log {log} \
            --thread {threads} \
        """

rule add_annotations:
    """Adds FORMAT annotations that were removed during processing with SHAPEIT4."""
    input: 
        vcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
        phased = config["results"] + "haplotypes/SHAPEIT4/{dataset}.{mode}.chr{chr}.vcf.gz",
    output:
        config["results"] + "haplotypes/SHAPEIT4/annotated/{dataset}.{mode}.chr{chr}.vcf.gz",
    threads: 1 
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools annotate {input.vcf} \
            -a {input.phased} \
            -c FORMAT/GT \
            -o {output} \
            -Oz \
        """

