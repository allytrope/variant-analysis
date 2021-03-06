"""Contain rules for determining relatedness among individuals.
These rules will generally be used after those in `variant_calling.smk`.
"""

CONFIG = config["relations"]

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
    input: vcf = CONFIG["vcf"],
    output: roh_pickle = config["results"] + "relatedness/roh/roh_poisson.pickle",  # Stores a pandas datafrane
            froh_pickle = config["results"] + "relatedness/roh/froh_poisson.pickle",  # Stores a dictionary
    conda: "../envs/scikit.yaml"
    script: "../scripts/runs_of_homozygosity.py"

# Under development
# rule scikit_allel_PCA:
#     """Perform principal component analysis of genotype data."""
#     input: vcf = CONFIG["vcf"],
#     output:
#     conda: "../envs/scikit.yaml"
#     script: "../scripts/PCA.py"
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

rule bed_to_ped:
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

