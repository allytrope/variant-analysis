"""Rules related to calculating admixture."""


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
    """Admixture calculated using RFMix. Expects phased data."""
    input:
        # BCFs are recommended
        query_VCF = config["results"] + "genotypes/pass/{dataset}.{subset}.SNP.chr{chr}.bcf",
        query_VCF_idx = config["results"] + "genotypes/pass/{dataset}.{subset}.SNP.chr{chr}.bcf.csi",
        ref_VCF = config["resources"] + "ref_vcf/12_Indian_12_Chinese.phased.autosomal.bcf",
        ref_VCF_idx = config["resources"] + "ref_vcf/12_Indian_12_Chinese.phased.autosomal.bcf.csi",
        sample_map = config["resources"] + "ref_vcf/sample_ancestry.tsv",
        genetic_map = config["resources"] + "genetic_map/autosomal.cM.genetic_map",
    params:
        output_base = config["results"] + "admixture/{dataset}.{subset}.RFMix.chr{chr}",
    output:
        # msp = lambda wildcards, params: params.output_base + "msp.tsv",
        # fb = lambda wildcards, params: params.output_base + "fb.tsv",
        msp = config["results"] + "admixture/{dataset}.{subset}.RFMix.chr{chr}.msp.tsv",
        fb = config["results"] + "admixture/{dataset}.{subset}.RFMix.chr{chr}.fb.tsv",
        Q = config["results"] + "admixture/{dataset}.{subset}.RFMix.chr{chr}.rfmix.Q",
        sis = config["results"] + "admixture/{dataset}.{subset}.RFMix.chr{chr}.sis.tsv",
    log: config["results"] + "admixture/{dataset}.{subset}.RFMix.chr{chr}.log",
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
        2> {log} \
        """

rule concat_msp:
    """Combine the .msp files for each chromosome of RFMix."""
    input:
        lambda wildcards: expand(config["results"] + "admixture/{dataset}.{subset}.RFMix.chr{chr}.msp.tsv",
            chr=[chrom for chrom in CHROMOSOMES if chrom not in ["X", "Y", "MT"]],
            dataset=wildcards.dataset,
            subset=wildcards.subset),
    output:
        concat = config["results"] + "admixture/{dataset}.{subset}.RFMix.all.msp.tsv",
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

# rule merge_Q:
#     input:
#         expand(config["results"] + "admixture/RFMix.chr{chr}.rfmix.Q", chr=[chrom for chrom in CHROMOSOMES if chrom not in ["X", "Y", "MT"]]),
#     output:
#     threads: 1
#     resources: nodes = 1
#     conda: "../envs/admixture.yaml"



    # shell: """
    #     TOTAL_LENGTH = grep -v "^#" {input.chr_stats} | grep "total-length" | grep "assembled-molecule" | cut -f 6 | awk '{{s+=$1}} END {{print s}}';
    #     for CHR in $CHROMOSOMES; do
    #         CHR_LENGTH = grep -v "^#" {input.chr_stats} | grep "total-length" | grep "assembled-molecule" |  awk 'BEGIN {FS="\t"} $2 == 1' | cut -f 6;
    #         PERCENT = $(echo "scale=5; $CHR_LENGTH / $TOTAL_LENGTH" | bc)
    #     done
    #     """
    # run:
    #     total_length = shell("grep -v '^#' {input.chr_stats} | grep 'total-length' | grep 'assembled-molecule' | cut -f 6 | awk '{{s+=$1}} END {{print s}}'")
    #     for file in input.Q:


## Unsuperived admixture using ADMIXTURE

# rule create_plink_files:
#    """Create PLINK .ped and .map files."""
#     input: vcf = CONFIG["vcf"],
#     output: ped = config["results"] + "plink/{dataset}.initial.ped",
#             map = config["results"] + "plink/{dataset}.map",
#     params: out = lambda wildcards: config["results"] + "plink/{dataset}".format(dataset=wildcards.dataset),
#     shell: "vcftools \
#                --vcf {input.vcf} \
#                --plink \
#                --out {params.out}; \
#             mv {params.out}.ped {params.out}.initial.ped"



# rule recode_plink_files:
#     """Convert base readings from {A,C,T,G} to {1,2}.
#     Not preferable. Ideally create from VCF."""
#     input:
#         lambda wildcards: expand("{prefix}.{ext}",
#             prefix=wildcards.prefix,
#             ext=["ped", "map"],
#         ),
#     output:
#         lambda wildcards: expand("{prefix}.recode12.{ext}",
#             prefix=wildcards.prefix,
#             ext=["ped", "map"],
#         ),
#     shell: """
#         plink \
#             --file {wildcards.prefix} \
#             --out {wildcards.prefix}.recode12 \
#             --recode12 \
#         """

# rule bcf_to_plink:
#     """Create .bed and associated PLINK files."""
#     input:
#         bcf = "{prefix}.bcf",
#     output:
#         lambda wildcards: expand("{prefix}.recode12.{ext}",
#             prefix=wildcards.prefix,
#             ext=["ped", "map"],
#         ),
#     shell: """
#         plink \
#             --file {input.bcf} \
#             --out {wildcards.prefix} \
#             --make-bed \
#             --recode12 \
#         """
#     shell: "plink \
#                 --bcf {input.vcf} \
#                 --recode12 \
#                 --make-bed \
#                 --out {config[results]}admixture/supervised/plink/{wildcards.sample}"

rule unsupervised_ADMIXTURE:
    """Cluster samples by ancestry in an unsupervised manner."""
    input:
        bed = config["results"] + "genotypes/pruned/plink/{dataset}.{subset}.{mode}.autosomal.bed",
        bim = config["results"] + "genotypes/pruned/plink/{dataset}.{subset}.{mode}.autosomal.bim",
        fam = config["results"] + "genotypes/pruned/plink/{dataset}.{subset}.{mode}.autosomal.fam",
    output:
        q = config["results"] + "admixture/ADMIXTURE/unsupervised/{dataset}.{subset}.{mode}.autosomal.{clusters}.Q",
        p = config["results"] + "admixture/ADMIXTURE/unsupervised/{dataset}.{subset}.{mode}.autosomal.{clusters}.P",
        log = config["results"] + "admixture/ADMIXTURE/unsupervised/{dataset}.{subset}.{mode}.autosomal.{clusters}.log",
    threads: 24
    conda: "../envs/admixture.yaml"
    # --cv | tee {output.out}; \
    shell: """
        admixture {input.bed} {wildcards.clusters} \
            -j{threads} \
            -s 899 \
            --cv > {output.log}; \
        mv {wildcards.dataset}.{wildcards.subset}.{wildcards.mode}.autosomal.{wildcards.clusters}.Q {output.q}; \
        mv {wildcards.dataset}.{wildcards.subset}.{wildcards.mode}.autosomal.{wildcards.clusters}.P {output.p}; \
        """


## Supervised ADMIXTURE
rule subset_ref:
    """Subset ref VCF."""
    input:
        ref_bcf = config["resources"] + "ref_vcf/Indian-Chinese.vcf.gz",
        ref_csi = config["resources"] + "ref_vcf/Indian-Chinese.vcf.gz.tbi",
    output:
        ref_bcf = config["resources"] + "ref_vcf/filtered/Indian-Chinese.chr{chr}.bcf",
    conda: "../envs/bio.yaml"
    shell: """
        bcftools view {input.ref_bcf} \
            -e 'F_MISSING>0.6' \
            --regions {wildcards.chr} \
        | bcftools annotate \
            -x INFO \
            -Ob \
            -o {output.ref_bcf} \
        """

rule intersect_ref:
    """Subset ref VCF to contain only SNVs in query VCF."""
    input:
        query_bcf = config["results"] + "genotypes/pass/{dataset}.{subset}.SNP.chr{chr}.bcf",
        query_csi = config["results"] + "genotypes/pass/{dataset}.{subset}.SNP.chr{chr}.bcf.csi",
        ref_bcf = config["resources"] + "ref_vcf/filtered/Indian-Chinese.chr{chr}.bcf",
        ref_csi = config["resources"] + "ref_vcf/filtered/Indian-Chinese.chr{chr}.bcf.csi",
    output:
        # query_bcf = config["results"] + "genotypes/pass/{dataset}.{subset}_Indian-Chinese_subset.SNP.chr{chr}.bcf",
        # ref_bcf = config["resources"] + "ref_vcf/Indian-Chinese_{dataset}{subset}_subset.chr{chr}.bcf",
        readme = config["results"] + "genotypes/isec/{dataset}.{subset}_Indian-Chinese.chr{chr}/README.txt",
        bcf0 = config["results"] + "genotypes/isec/{dataset}.{subset}_Indian-Chinese.chr{chr}/0000.bcf",
        bcf1 = config["results"] + "genotypes/isec/{dataset}.{subset}_Indian-Chinese.chr{chr}/0001.bcf",
        sites = config["results"] + "genotypes/isec/{dataset}.{subset}_Indian-Chinese.chr{chr}/sites.txt",
    params:
        prefix = lambda wildcards, output: "/".join(output.readme.split("/")[0:-1]),
    threads: 2
    resources:
        nodes = 2
    conda: "../envs/bio.yaml"
    shell: """
        bcftools isec \
            {input.query_bcf} \
            {input.ref_bcf} \
            --collapse none \
            -n=2 \
            -Ob \
            -p {params.prefix} \
        """

rule merge_ref:
    """Create reference VCF (to be turned into PLINK BED) for supervised ADMIXTURE.
    This must have both the samples to be queried and the reference samples."""
    input:
        bcf0 = config["results"] + "genotypes/isec/{dataset}.{subset}_Indian-Chinese.chr{chr}/0000.bcf",
        bcf1 = config["results"] + "genotypes/isec/{dataset}.{subset}_Indian-Chinese.chr{chr}/0001.bcf",
    output:
        bcf = config["results"] + "genotypes/pass/{dataset}.{subset}_Indian-Chinese_merged.SNP.chr{chr}.bcf",
    conda: "../envs/bio.yaml"
    shell: """
        bcftools merge {input.bcf0} {input.bcf1} \
            -Ob \
            -o {output.bcf} \
        """

rule make_pop_ADMIXTURE:
    """Create .pop file with same order as individuals in .fam.
    This is for ADMIXTURE. This file has no header either.
    Each line is an origin."""
    input:
        #fam = "{prefix}/{name}.fam",
        #fam = "{prefix}/{dataset}.{subset}_Indian-Chinese_merged.{mode}.autosomal.fam",
        fam = "{prefix}/{dataset}.{subset}.{mode}.autosomal.fam",
        origins = config["resources"] + "pop/founder_origins.tsv",  #$ For founders as ancestral populations
        #origins = config["resources"] + "ref_vcf/Indian-Chinese_ancestry.tsv",  # For Indian and Chinese as ancestral populatiosn
    output:
        #population = "{prefix}/{name}.pop",
        #population = "{prefix}/{dataset}.{subset}_Indian-Chinese_merged.{mode}.autosomal.pop",
        population = "{prefix}/{dataset}.{subset}.{mode}.autosomal.pop",
    conda: "../envs/bio.yaml"
    #-f 'Indiv;Id' \
    shell: """
        csvtk join \
            <(cat <(echo -e 'Fam\tIndiv\tSire\tDam\tSex\tPhenotype') <(sed 's/ /\t/g' {input.fam})) \
            {input.origins} \
            -t \
            -f 'Indiv' \
            --left-join \
        | csvtk cut \
            -t \
            -f Indiv,Origin \
        | csvtk cut \
            -t \
            -f Origin \
        | sed '1d' \
        > {output.population} \
        """

rule supervised_ADMIXTURE:
    """Cluster samples by ancestry in supervised manner.
    Requires .pop file and related PLINK files in appropriate directory as below.
    These are not generated by any rules and should have truth samples, those with known ancestries."""
    # Make {subset} end with `_Indian-Chinese_merged` in order to trigger above rules.
    input:
        bed = config["results"] + "genotypes/pruned/plink/{dataset}.{subset}.{mode}.autosomal.bed",
        bim = config["results"] + "genotypes/pruned/plink/{dataset}.{subset}.{mode}.autosomal.bim",
        fam = config["results"] + "genotypes/pruned/plink/{dataset}.{subset}.{mode}.autosomal.fam",
        population = config["results"] + "genotypes/pruned/plink/{dataset}.{subset}.{mode}.autosomal.pop",
    # params:
    #     prefix = lambda output: ".".join(output.q.split(".")[0:-1]),
    output:
        q = config["results"] + "admixture/ADMIXTURE/supervised/{dataset}.{subset}.{mode}.autosomal.{clusters}.Q",
        p = config["results"] + "admixture/ADMIXTURE/supervised/{dataset}.{subset}.{mode}.autosomal.{clusters}.P",
    log: config["results"] + "admixture/ADMIXTURE/supervised/{dataset}.{subset}.{mode}.autosomal.{clusters}.log",
    threads: 24
    conda: "../envs/admixture.yaml"
    # Do I just remove `{wildcards.clusters}` from being the second positional argument?
    shell: """
        admixture {input.bed} {wildcards.clusters} \
            --supervised \
            -j{threads} \
            -s 899 \
            --cv > {log}; \
        mv {wildcards.dataset}.{wildcards.subset}.{wildcards.mode}.autosomal.{wildcards.clusters}.Q {output.q}; \
        mv {wildcards.dataset}.{wildcards.subset}.{wildcards.mode}.autosomal.{wildcards.clusters}.P {output.p}; \
        """


# rule create_unsupervised_ancestry_tsv:
#     """Remove added knowns from .Q and add sample IDs. This makes a more human-readable version of .Q."""
#     input: ped = config["results"] + "plink/recode12/{dataset}_recode12.ped",
#            q = config["results"] + "admixture/unsupervised/{dataset}.{clusters}.Q",
#     output: ancestries = config["results"] + "admixture/unsupervised/{dataset}.{clusters}.admixture.tsv",
#     # `echo` adds header line
#     # `sed` replaces spaces with tabs as delimiter.
#     # `cut` takes first second column of .ped
#     # `paste` combines the outputs of `sed` and `cut`
#     # `head` removes known truth lines from end of file
#     # GENERALIZE HEADER CREATION in `echo` command
#     # GENERALIZE REMOVAL OF KNOWN TRUTH LINES in `head` command
#     shell: "echo -e 'sample_id\tCluster_1\tCluster_2' > {output}; \
#             TMP=delimited.tmp; \
#             sed 's/ /\t/g' {input.q} > $TMP; \
#             cut {input.ped} -d ' ' -f 2 \
#             | paste - $TMP >> {output}; \
#             rm $TMP;"