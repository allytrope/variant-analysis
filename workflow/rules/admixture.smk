## Must be modified for different runs

rule convert_to_cM:
    """Convert cM/Mb to cM."""
    input:
        recomb_rates_chrom = "/master/abagwell/workspace/linkage_disequalibrium/chr{chr}.bg",
    output:
        genetic_map = "/master/abagwell/workspace/linkage_disequalibrium/chr{chr}.cM.bg",
    #START=$(grep -P -m 1 'chr{wildcards.chr}\t' chr1.bg | cut -f 2); \
    shell: """
        awk -v START=0 'BEGIN {{PREV_START=START}} {{$4=(($2-PREV_START)*$4*0.000001)+PREV_CM; printf "%s %s %s %0.9f\n", $1, $2, $3, $4; PREV_START=$2; PREV_CM=$4}}' {input.recomb_rates_chrom} \
        > {output.genetic_map}
        """

rule liftover_recomb_rates:
    """Liftover BEDGRAPH by chromosome."""
    input:
        recomb_rates = "/master/abagwell/workspace/linkage_disequalibrium/chr{chr}.cM.bg",
        chain = "/master/abagwell/workspace/linkage_disequalibrium/rheMac8ToRheMac10.over.chain.gz",
    output:
        remapped = "/master/abagwell/workspace/linkage_disequalibrium/chr{chr}.cM.remapped.bg",
        unmapped = "/master/abagwell/workspace/linkage_disequalibrium/chr{chr}.cM.unmapped.bg",
    shell: """
        liftOver {input.recomb_rates} {input.chain} {output.remapped} {output.unmapped}; \
        """

rule RFMix:
    """Admixture calculated using RFMix."""
    input:
        # BCFs are recommended
        #query_VCF = config["results"] + "haplotypes/SHAPEIT4/with_scaffold/all_samples.SNP.chr{chr}.vcf.gz",
        query_VCF = config["results"] + "haplotypes/SHAPEIT5/all_samples.SNP.chr{chr}.vcf.gz",
        ref_VCF = config["resources"] + "ref_vcf/12_Indian_12_Chinese.phased.chr{chr}.vcf.gz",
        ref_VCF_idx = config["resources"] + "ref_vcf/12_Indian_12_Chinese.phased.chr{chr}.vcf.gz.tbi",
        sample_map = config["resources"] + "ref_vcf/sample_ancestry.tsv",
        genetic_map =config["resources"] + "genetic_map/chr{chr}.cM.genetic_map",
    params:
        output_base = config["results"] + "admixture/RFMix.chr{chr}",
    output:
        # msp = lambda wildcards, params: params.output_base + "msp.tsv",
        # fb = lambda wildcards, params: params.output_base + "fb.tsv",
        msp = config["results"] + "admixture/RFMix.chr{chr}.msp.tsv",
        fb = config["results"] + "admixture/RFMix.chr{chr}.fb.tsv",
        Q = config["results"] + "admixture/RFMix.chr{chr}.rfmix.Q",
        sis = config["results"] + "admixture/RFMix.chr{chr}.sis.tsv",
    threads: 1
    resources: nodes = 1
    conda: "../envs/admixture.yaml"
    shell: """
        rfmix \
            -f {input.query_VCF} \
            -r {input.ref_VCF} \
            -m {input.sample_map} \
            -g {input.genetic_map} \
            -o {params.output_base} \
            --chromosome={wildcards.chr} \
            -e 3 \
            --reanalyze-reference \
        """

rule concat_msp:
    """Combine the .msp files for each chromosome of RFMix."""
    input:
        expand(config["results"] + "admixture/RFMix.chr{chr}.msp.tsv", chr=[chrom for chrom in CHROMOSOMES if chrom not in ["X", "Y", "MT"]]),
    output:
        concat = config["results"] + "admixture/RFMix.all.msp.tsv",
    threads: 1
    resources: nodes = 1
    conda: "../envs/admixture.yaml"
    shell: """
        FILES=({input}); \
        head -n 2 ${{FILES[0]}} > {output.concat}; \
        for FILE in "${{FILES[@]}}"; do \
            sed '1d;2d' $FILE >> {output.concat}; \
        done \
        """