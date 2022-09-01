"""Contain rules for determining relatedness among individuals.
These rules will generally be used after those in `variant_calling.smk`.
"""

CONFIG = config["relations"]

## Calculate relatedness with KING


rule concat_chromosomes_skip_SHAPEIT4:
    """Concatenate chromosomes."""
    input:
        vcfs = lambda wildcards: expand(config["results"] + "genotypes/filtered/{dataset}.SNP.chr{chr}.filtered.vcf.gz",  # Taking hard filtered .vcf,
            dataset=wildcards.dataset,
            chr=[i for i in range(1, 21)] + ['X']),
    output:
        config["results"] + "genotypes/filtered/{dataset}.SNP.filtered.vcf.gz"
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools concat {input.vcfs} \
            -o {output} \
            -Oz \
        """

rule concat_chromosomes:
    """Concatenate chromosomes."""
    input:
        vcfs = lambda wildcards: expand(config["results"] + "haplotypes/SHAPEIT4/{dataset}.SNP.chr{chr}.phased.vcf.gz",
            dataset=wildcards.dataset,
            chr=[i for i in range(1, 21)] + ['X']),
    output:
        config["results"] + "haplotypes/SHAPEIT4/{dataset}.SNP.phased.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        bcftools concat {input.vcfs} \
            -o {output} \
            -Oz \
        """

rule create_combined_binary_plink_file:
    """Create .bed and associated PLINK files."""
    input: 
        #vcf = config["results"] + "haplotypes/SHAPEIT4/{dataset}.SNP.phased.vcf.gz",  # Contains all chromosomes
        vcf = config["results"] + "genotypes/filtered/{dataset}.SNP.filtered.min-ac_3.vcf.gz",
        parents_table = config["parents_table"],
        sex_table = config["sex_table"],
    #output: expand("{results}admixture/supervised/plink/{sample}.{ext}", results=config["results"], sample="{wildcards.sample}", ext=["bed", "bim", "fam"])
    output:
        bed = config["results"] + "kinship/plink/{dataset}.bed",
        bim = config["results"] + "kinship/plink/{dataset}.bim",
        fam = config["results"] + "kinship/plink/{dataset}.fam",
    params:
        out_path = config["results"] + "kinship/plink",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        plink \
            --vcf {input.vcf} \
            --update-parents {input.parents_table} \
            --update-sex {input.sex_table} \
            --recode12 \
            --make-bed \
            --out {params.out_path}/{wildcards.dataset} \
        """

rule set_as_same_family_plink:
    """Set the FID field for all individuals to 1."""
    input:
        bed = config["results"] + "kinship/plink/{dataset}.bed",
        bim = config["results"] + "kinship/plink/{dataset}.bim",
        fam = config["results"] + "kinship/plink/{dataset}.fam",
        update_fids = config["update_fids"],
    params: path = config["results"] + "kinship/plink",
    output:
        bed = config["results"] + "kinship/plink/{dataset}.same_fid.bed",
        bim = config["results"] + "kinship/plink/{dataset}.same_fid.bim",
        fam = config["results"] + "kinship/plink/{dataset}.same_fid.fam",
    shell: """
        plink \
            --bfile {params.path}/{wildcards.dataset} \
            --update-ids {input.update_fids} \
            --make-bed \
            --out {params.path}/{wildcards.dataset}.same_fid \
        """

rule estimate_relatedness_king:
    """Estimate relatedness between individuals."""
    input: 
        bed = config["results"] + "kinship/plink/GBS_WES_WGS.same_fid.bed",
        bim = config["results"] + "kinship/plink/GBS_WES_WGS.same_fid.bim",
        fam = config["results"] + "kinship/plink/GBS_WES_WGS.same_fid.fam",
    output:
        kinship = config["results"] + "kinship/KING/king.kin",
    params:
        prefix = config["results"] + "kinship/KING/king",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        king \
            -b {input.bed} \
            --kinship \
            --prefix {params.prefix} \
        """


## Per chromosome

rule create_combined_binary_plink_file_per_chr:
    """Create .bed and associated PLINK files."""
    input:
        vcf = config["results"] + "genotypes/filtered/{dataset}.SNP.chr{chr}.filtered.vcf.gz",
        #vcf = config["results"] + "haplotypes/SHAPEIT4/{dataset}.SNP.chr{chr}.phased.vcf.gz",  # Contains all chromosomes
        parents_table = config["parents_table"],
        sex_table = config["sex_table"],
    #output: expand("{results}admixture/supervised/plink/{sample}.{ext}", results=config["results"], sample="{wildcards.sample}", ext=["bed", "bim", "fam"])
    output:
        bed = config["results"] + "kinship/plink/{dataset}.chr{chr}.bed",
        bim = config["results"] + "kinship/plink/{dataset}.chr{chr}.bim",
        fam = config["results"] + "kinship/plink/{dataset}.chr{chr}.fam",
    params:
        out_path = config["results"] + "kinship/plink",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        plink \
            --vcf {input.vcf} \
            --update-parents {input.parents_table} \
            --update-sex {input.sex_table} \
            --recode12 \
            --make-bed \
            --out {params.out_path}/{wildcards.dataset}.chr{wildcards.chr} \
        """

rule set_as_same_family_plink_per_chr:
    """Set the FID field for all individuals to 1."""
    input:
        bed = config["results"] + "kinship/plink/{dataset}.chr{chr}.bed",
        bim = config["results"] + "kinship/plink/{dataset}.chr{chr}.bim",
        fam = config["results"] + "kinship/plink/{dataset}.chr{chr}.fam",
        update_fids = config["update_fids"],
    params: path = config["results"] + "kinship/plink",
    output:
        bed = config["results"] + "kinship/plink/{dataset}.chr{chr}.same_fid.bed",
        bim = config["results"] + "kinship/plink/{dataset}.chr{chr}.same_fid.bim",
        fam = config["results"] + "kinship/plink/{dataset}.chr{chr}.same_fid.fam",
    shell: """
        plink \
            --bfile {params.path}/{wildcards.dataset}.chr{wildcards.chr} \
            --update-ids {input.update_fids} \
            --make-bed \
            --out {params.path}/{wildcards.dataset}.chr{wildcards.chr}.same_fid \
        """

rule estimate_relatedness_king_per_chr:
    """Estimate relatedness between individuals."""
    input: 
        bed = config["results"] + "kinship/plink/GBS_WES_WGS.chr{chr}.same_fid.bed",
        bim = config["results"] + "kinship/plink/GBS_WES_WGS.chr{chr}.same_fid.bim",
        fam = config["results"] + "kinship/plink/GBS_WES_WGS.chr{chr}.same_fid.fam",
    output:
        kinship = config["results"] + "kinship/KING/king.chr{chr}.kin",
    params:
        prefix = config["results"] + "kinship/KING/king.chr{chr}",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        king \
            -b {input.bed} \
            --kinship \
            --prefix {params.prefix} \
        """




## Inbreeding

rule inbreeding_coefficient:
    """Calculate inbreeding coefficient, which is a measure of heterozygosity per-individual."""
    input: vcf = config["results"] + "haplotypes/SHAPEIT4/GBS_WES_WGS.SNP.phased.vcf.gz",
    output: het = config["results"] + "kinship/het/GBS_WES_WGS.het",
    #params: het = config["results"] + "kinship/het/GBS_WES_WGS",
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        vcftools \
            --gzvcf {input.vcf} \
            --het \
            --stdout \
        > {output.het} \
        """



rule estimate_relatedness_lcmlkin:
    """Estimate relatedness between individuals using maximum likelihood."""
    input: vcf = CONFIG["vcf"],
           founders = CONFIG["founders"],
    output: config["results"] + "relatedness/recalibrated_joint_call.relate",
    threads: 8  # max of 8
    shell: "lcmlkin \
                -i {input.vcf} \
                -o {output} \
                -g all \
                -l phred \
                -u {input.founders} \
                -t {threads}"

# Find Fst using scikit-allel
rule scikit_allel_fst:
    """Calculate pairwise Fst between populations within non-overlapping windows.
    Output file stores data from a pandas dataframe to be easily graphed with a Python graphing library like seaborn."""
    input: vcf = CONFIG["vcf"],
    output: #config["results"] + "relatedness/fst/created_scikit_fst.txt",
            #hdf5 = config["results"] + "relatedness/fst/fsts.hdf5",
            #csv = config["results"] + "relatedness/fst/fsts.csv",
            pickle = config["results"] + "relatedness/fst/fsts.pickle",
    conda: "../envs/scikit.yaml"
    params: subpops = CONFIG["fst"]["pops"],
    script: "../scripts/fst.py"

rule scikit_allel_diversity:
    """Estimate nucleotide diversity in windows. Output files are formatted to be graphed with the R tool chromoMap."""
    input: vcf = CONFIG["vcf"],
    output: annotations = config["results"] + "relatedness/diversity/diversity.tsv",
            chromosomes = config["results"] + "relatedness/diversity/chromosomes.tsv",
    conda: "../envs/scikit.yaml"
    script: "../scripts/diversity.py"

rule scikit_allel_ROH:
    """Calculate runs of homozygosity."""
    input:
        vcf = config["results"] + "haplotypes/SHAPEIT4/{dataset}.{mode}.chr{chr}.phased.vcf.gz",
    #vcf = CONFIG["vcf"],
    output:
        roh_pickle = config["results"] + "relatedness/roh/{dataset}.{mode}.chr{chr}.roh_poisson.pickle",  # Stores a pandas datafrane
        froh_pickle = config["results"] + "relatedness/roh/{dataset}.{mode}.chr{chr}.froh_poisson.pickle",  # Stores a dictionary
    threads: 1
    resources: nodes = 1
    conda: "../envs/scikit.yaml"
    script: "../scripts/runs_of_homozygosity.py"

# Under development
# rule scikit_allel_PCA:
#     """Perform principal component analysis of genotype data."""
#     input: vcf = CONFIG["vcf"],
#     output:
#     conda: "../envs/scikit.yaml"
#     script: "../scripts/pca.py"
# --------------
# Unsupervised admixture


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


'''
rule recode_plink_files:
    """Convert base readings from {A,C,T,G} to {1,2}."""
    input: ped=config["results"] + "plink/{dataset}.ped",
           map=config["results"] + "plink/{dataset}.map"
    output: ped=config["results"] + "plink/recode12/{dataset}_recode12.ped",
            map=config["results"] + "plink/recode12/{dataset}_recode12.map"
    shell: "plink \
                --file  " + config["results"] + " plink/{wildcards.dataset} \
                --out " + config["results"] + " plink/recode12/{wildcards.dataset}_recode12 \
                --recode12"
'''

rule admixture:
    """Cluster samples by ancestry in an unsupervised manner."""
    input: ped = config["results"] + "plink/recode12/{dataset}_recode12.ped",
           map = config["results"] + "plink/recode12/{dataset}_recode12.map",
    output: q = config["results"] + "admixture/unsupervised/{dataset}.{clusters}.Q",
            p = config["results"] + "admixture/unsupervised/{dataset}.{clusters}.P",
            out = config["results"] + "admixture/unsupervised/{dataset}.{clusters}.out",
    threads: 24
    conda: "../envs/bio.yaml"
    shell: "admixture {input.ped} {wildcards.clusters} \
                -j{threads} \
                --cv | tee {output.out}; \
            mv {wildcards.dataset}_recode12.{wildcards.clusters}.Q {output.q}; \
            mv {wildcards.dataset}_recode12.{wildcards.clusters}.P {output.p}"

rule create_unsupervised_ancestry_tsv:
    """Remove added knowns from .Q and add sample IDs. This makes a more human-readable version of .Q."""
    input: ped = config["results"] + "plink/recode12/{dataset}_recode12.ped",
           q = config["results"] + "admixture/unsupervised/{dataset}.{clusters}.Q",
    output: ancestries = config["results"] + "admixture/unsupervised/{dataset}.{clusters}.admixture.tsv",
    # `echo` adds header line
    # `sed` replaces spaces with tabs as delimiter.
    # `cut` takes first second column of .ped
    # `paste` combines the outputs of `sed` and `cut`
    # `head` removes known truth lines from end of file
    # GENERALIZE HEADER CREATION in `echo` command
    # GENERALIZE REMOVAL OF KNOWN TRUTH LINES in `head` command
    shell: "echo -e 'sample_id\tCluster_1\tCluster_2' > {output}; \
            TMP=delimited.tmp; \
            sed 's/ /\t/g' {input.q} > $TMP; \
            cut {input.ped} -d ' ' -f 2 \
            | paste - $TMP >> {output}; \
            rm $TMP;"

# ---------------
# Supervised admixture
#rule merge_vcfs:
#	"""Combine vcfs."""
#	input:
#	output:

#rule filter_vcf:
#	input:
#	output:

rule create_binary_plink_files:
    """Create .bed and associated PLINK files."""
    input: vcf = config["results"] + "admixture/supervised/{sample}.vcf.gz",
    #output: expand("{results}admixture/supervised/plink/{sample}.{ext}", results=config["results"], sample="{wildcards.sample}", ext=["bed", "bim", "fam"])
    output: bed = config["results"] + "admixture/supervised/plink/{sample}.bed",
            bim = config["results"] + "admixture/supervised/plink/{sample}.bim",
            fam = config["results"] + "admixture/supervised/plink/{sample}.fam",
    conda: "../envs/bio.yaml"
    shell: "plink \
                --vcf {input.vcf} \
                --recode12 \
                --make-bed \
                --out {config[results]}admixture/supervised/plink/{wildcards.sample}"

rule supervised_admixture:
    """Cluster samples by ancestry in supervised manner.
    Requires .pop file and related PLINK files in appropriate directory as below.
    These are not generated by any rules and should have truth samples, those with known ancestries."""
    input: bed = config["results"] + "admixture/supervised/plink/{dataset}.bed",
           bim = config["results"] + "admixture/supervised/plink/{dataset}.bim",
           fam = config["results"] + "admixture/supervised/plink/{dataset}.fam",
           population = config["results"] + "admixture/supervised/plink/{dataset}.pop",
    output: q = config["results"] + "admixture/supervised/{dataset}.{clusters}.Q",
            p = config["results"] + "admixture/supervised/{dataset}.{clusters}.P",
    threads: 24
    conda: "../envs/bio.yaml"
    shell: "admixture {input.bed} {wildcards.clusters} \
                --supervised \
                -j{threads}; \
            mv {wildcards.dataset}.{wildcards.clusters}.Q {output.q}; \
            mv {wildcards.dataset}.{wildcards.clusters}.P {output.p}"

rule bed_to_ped_whole_genome:
    """Convert PLINK's binary .bed into nonbinary .ped."""
    input: bed = config["results"] + "admixture/supervised/plink/{dataset}.bed",
           bim = config["results"] + "admixture/supervised/plink/{dataset}.bim",
           fam = config["results"] + "admixture/supervised/plink/{dataset}.fam",
    output: 
           ped = config["results"] + "admixture/supervised/plink/{dataset}.ped",
           map = config["results"] + "admixture/supervised/plink/{dataset}.map",
    conda: "../envs/bio.yaml"
    shell: "plink \
                --bfile " + config["results"] + "admixture/supervised/plink/{wildcards.dataset} \
                --recode \
                --out " + config["results"] + "admixture/supervised/plink/{wildcards.dataset}"

rule create_supervised_ancestry_tsv:
    """Remove added knowns from .Q and add sample IDs. This makes a more human-readable version of .Q."""
    input: ped = config["results"] + "admixture/supervised/plink/{dataset}.ped",
           q = config["results"] + "admixture/supervised/{dataset}.{clusters}.Q",
           population = config["results"] + "admixture/supervised/plink/{dataset}.pop",
    output: ancestries = config["results"] + "admixture/supervised/{dataset}.{clusters}.admixture.tsv",
    # `echo` adds header line
    # `sed` replaces spaces with tabs as delimiter.
    # `cut` takes first second column of .ped
    # `paste` combines the outputs of `sed` and `cut`
    # `head` removes known truth lines from end of file

    # GENERALIZE HEADER CREATION in `echo` command
    # GENERALIZE REMOVAL OF KNOWN TRUTH LINES in `head` command
    shell: "echo -e 'Sample_id\tChinese\tIndian' > {output}; \
            TMP=delimited.tmp; \
            TMP2=pasted.tmp; \
            sed 's/ /\t/g' {input.q} > $TMP; \
            cut {input.ped} -d ' ' -f 2 \
            | paste - $TMP >> $TMP2; \
            tail -n 674 $TMP2 >> {output}; \
            rm $TMP; rm $TMP2"

# Requires some modification depending on data.
rule graph_admixture:
    """Graph ancestry of samples."""
    input: config["results"] + "admixture/supervised/{dataset}.{clusters}.admixture.tsv",
    output: config["results"] + "admixture/supervised/{dataset}.{clusters}.Q.png",
    run:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns

        data = pd.read_table(str(input))
        sns.set_theme()
        ax = sns.histplot(data=data, bins=20, multiple="stack")
        ax.set(title="Ancestral Admixture in SNPRC Rhesus Macaques", xlabel="Ancestry ratio", xlim=(0.5, 1), ylabel="Number of individuals")
        plt.savefig(str(output))


rule find_AIMs:
    """Find ancestry informative markers, variants that can be used to predict ancestry. Uses a threshold value.
    If say threshold is 0.85, that means, it requires a minimum of 0.85 frequency in one cluster
    and maximum of 0.15 in others to be considered an AIM. Currently works only when 2 clusters are involved."""
    input: p=config["results"] + "admixture/supervised/{dataset}.{clusters}.P",
           map=config["results"] + "admixture/supervised/plink/{dataset}.map",
    output: config["results"] + "aims/{dataset}.{clusters}.aims",
    params: max_diff = CONFIG["AIMs_max_diff"],  # Frequency below which to not consider as an AIM.
    # Select chromosome and position columns from map and combine with .P by row.
    # Then filter for rows that fulfill the thresholds.
    shell: "awk '{{print $1,$4}}' {input.map} \
            | paste - {input.p} \
            | awk '$4 - $3 > {params.max_diff} || $3 - $4 > {params.max_diff} {{print $0}}' - \
            > {output}"

# rule create_chromosome_file:
#     """Chromosome file for `chromoMap`."""
#     input: "",
#     output: "",
#     shell: "awk 'BEGIN {{ID = 1}} {{print $1,1,$2,$2+1}}' {input} > {output}"

# rule create_annotation_file:
#     """Annotation file for `chromoMap`."""
#     input: config["results"] + "aims/{dataset}.{clusters}.aims",
#     output: config["results"] + "aims",
#     shell: "awk 'BEGIN {{ID = 1}} {{print \"SNP\"ID,$1,$2,$2+1; ID += 1}}' {input} > {output}"

