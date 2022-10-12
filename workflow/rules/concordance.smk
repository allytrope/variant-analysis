"""Rules related to determining concordance between samples."""


rule separate_sample:
    """Create single-sample VCF file."""
    input:
        #joint_vcf = config["results"] + "genotypes/filtered/{dataset}_labels.{mode}.chr{chr}.filtered.vcf.gz",
        joint_vcf = config["results"] + "/genotypes/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
    output:
        indiv_vcf = temp(config["results"] + "concordance/indiv_vcfs/{dataset}.{sample}.{mode}.chr{chr}.vcf.gz"),
    shell: """
        bcftools view {input.joint_vcf} \
            -s {wildcards.sample} \
            -Oz \
            -o {output.indiv_vcf} \
        """

rule genotype_concordance:
    """Calculate genotype concordance between two samples."""
    input:
        vcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
    output:
        concordance = config["results"] + "concordance/matches/{dataset}.{sample1}_{sample2}.{mode}.chr{chr}.concordance.tsv",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    # Reduce columns
    # Remove header
    # Keep only GT columns
    # Replace "|" with "/"
    # Swap any "1/0" with "0/1"
    # Count matches and mismatches
    shell: """
        bcftools annotate {input.vcf} \
            -x INFO,FORMAT \
        | bcftools view \
            -s {wildcards.sample1},{wildcards.sample2} \
        | bcftools view \
            -i 'F_MISSING=0' -H \
        | cut -f 10-11 \
        | sed 's/|/\//g' \
        | sed 's/1\/0/0\/1/g' \
        | awk -v OFS='\t' \
            'BEGIN {{print "SAMPLE1","SAMPLE2","MATCHES","MISMATCHES","CONCORDANCE"}} \
            {{if ($1 == $2) MATCH+=1; else MISMATCH+=1}} \
            END {{CONCORDANCE = MATCH/(MATCH + MISMATCH); print "{wildcards.sample1}","{wildcards.sample2}",MATCH,MISMATCH,CONCORDANCE}}' \
        > {output.concordance} \
        """

rule combine_concordances:
    """Combine concordances information from different samples."""
    input:
        lambda wildcards: expand(config["results"] + "concordance/matches/unfiltered/{dataset}.WGS{animal_id}_GBS{animal_id}.{mode}.chr{chr}.unfiltered.concordance.tsv",
            dataset=wildcards.dataset,
            animal_id=SAMPLE_NAMES,
            mode=wildcards.mode,
            chr=wildcards.chr),
    output: config["results"] + "concordance/{dataset}.{mode}.chr{chr}.unfiltered.concordance.tsv",
    # Take header
    # Take each file, remove header, and add 0's to empty field
    shell: """
            cat {input[0]} | head -n 1 > {output}; \
            for FILE in {input}; do \
                cat $FILE | sed '1d' | sed 's/\t\t/\t0\t/g' >> {output}; \
            done \
        """

rule count_non_missing:
    """Generate a table of counts."""
    input:
        vcf = config["results"] + "/genotypes/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
    output:
        called_count = config["results"] + "concordance/counts/{sample}.txt",
    shell: """
        bcftools view {input.vcf} \
            -s {wildcards.sample} \
        | bcftools view \
            -i "F_MISSING=0" -H \
        | wc -l \
        > {output.called_count} \
        """

rule tabulate_non_missing:
    """Generate a table of counts."""
    input:
        called_counts = expand(config["results"] + "concordance/counts/{sample}.txt", sample=SAMPLE_NAMES),
    output:
        combined_counts = config["results"] + "concordance/counts.list",
        tsv = config["results"] + "concordance/counts.tsv",
    shell: """
        cat {input.called_counts} > {output.combined_counts}; \
        paste {config[samples]} {output.combined_counts} > {output.tsv} \
        """