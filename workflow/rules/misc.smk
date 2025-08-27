# """Micellaneous rules."""
# rule bcf_to_vcfgz:
#     input: bcf = "{path}/{name}.bcf"
#     output: vcfgz = "{path}/{name}.vcf.gz"
#     conda: "../envs/common.yaml"
#     shell: """
#         bcftools view {input} \
#             -Oz \
#             -o {output.vcfgz} \
#         """

# rule vcfgz_to_bcf:
#     input: vcfgz = "{path}/{name}.vcf.gz"
#     output: bcf = "{path}/{name}.bcf"
#     conda: "../envs/common.yaml"
#     shell: """
#         bcftools view {input} \
#             -Ob \
#             -o {output.bcf} \
#         """

rule slivar_de_novo:
    input:
        bcf = config["results"] + "genotypes/pass/{dataset}.all.both.chr{chr}.bcf",
        ped = config["resources"] + "pedigree/plink.ped",
    output:
        bcf = config["results"] + "slivar_trios/de_novo/{dataset}.chr{chr}.bcf",
    params:
        slivar_functions = "/master/abagwell/tools/slivar-functions.js",
    threads: 1
    resources: nodes = 1
    conda:  "../envs/ensembl.yaml",
    shell: """
        slivar expr \
            --js {params.slivar_functions} \
            --vcf {input.bcf} \
            --ped {input.ped} \
            --trio "DENOVO:denovo(kid, dad, mom)" \
        | grep -E '^#|DENOVO' \
        | bcftools view -Ob -o {output.bcf} \
        """

rule vep_slivar_trio:
    """Annotate using ENSEMBL VEP."""
    wildcard_constraints:
        dam="[0-9]+"
    input:
        bcfs = lambda wildcards: expand(config["results"] + "slivar_trios/de_novo/{child}_{sire}_{dam}.chr{chr}.bcf",
            child=wildcards.child,
            sire=wildcards.sire,
            dam=wildcards.dam,
            chr=[str(num) for num in list(range(1,21))]),
        fasta = config["ref_fasta"],
        gtf = config["gtf3"],
        gtf_idx = config["gtf3"] + ".tbi",
    output:
        vcf = config["results"] + "slivar_trios/impact/{child}_{sire}_{dam}.vcf",
        stats = config["results"] + "slivar_trios/impact/{child}_{sire}_{dam}.html",
    params:
        species = "macaca_mulatta",
    threads: 1
    resources: nodes = 1
    conda:  "../envs/ensembl.yaml",
    #             --compress_output bgzip \

    # | bcftools view \
    # -i 'IMPACT="HIGH" | IMPACT="MODERATE"' \
    shell: """
        bcftools concat {input.bcfs} \
            -Ou \
        | bcftools view \
            -f PASS \
            -Ov \
        | vep \
            --fasta {input.fasta} \
            --format vcf \
            --gtf {input.gtf} \
            --output_file STDOUT \
            --stats_file {output.stats} \
            --species {params.species} \
            --vcf \
        | bcftools +split-vep \
            -c 'IMPACT' \
            -Ov \
            -o {output.vcf} \
        """

## List chromosomes or subsets

rule list_contigs:
    """Find all contigs from headers of reference genome"""
    input:
        ref = config["ref_fasta"],
    output:
        config["resources"] + "ref_fna/contigs.list",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        zcat {input.ref} | grep "^>" | cut -c 2- | cut -d " " -f 1 > {output}
        """

rule list_chromosomes:
    """Find all chromosomes.
    
     Each line of the output file is just one chromosome's name.
     The search only keeps numbered chromosomes (not those prefixed with "chr" or any other letters)
     as well as X, Y, and MT. Unplaced contigs are ignored."""
    input:
        contigs = config["resources"] + "ref_fna/contigs.list",
    output:
        config["resources"] + "ref_fna/chromosomes.list",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        cat {input.contigs} | grep -x -E "^[0-9]+|X|Y|MT" > {output}
        """

rule list_autosomes:
    """Find all autosomes. That is, only the numbered chromosomes. Does not include unplaced contigs."""
    input:
        chromosomes = config["resources"] + "ref_fna/chromosomes.list",
    output:
        config["resources"] + "ref_fna/autosomes.list",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        cat {input.chromosomes} | grep -vxE "X|Y|MT" > {output}
        """

################





rule concat_chromosomes:
    """Concatenate chromosomes."""
    input:
        vcfs = lambda wildcards: expand(config["results"] + "{path}/{dataset}.{mode}.chr{chr}.vcf.gz",
            path=wildcards.path,
            dataset=wildcards.dataset,
            mode=wildcards.mode,
            chr=CHROMOSOMES),
        chromosomes = config["resources"] + "ref_fna/chromosomes.list",
    output:
        vcf = config["results"] + "{path}/{dataset}.{mode}.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools concat {input.vcfs} \
            -o {output.vcf} \
            -Oz \
        """

rule bcfs_to_autosomal_bcf:
    """Concatenate autosomes."""
    wildcard_constraints:
        ext = "bcf",
    input:
        bcfs = lambda wildcards: expand(config["results"] + "{path}.chr{chr}.{ext}",
            path=wildcards.path,
            chr=[chrom for chrom in CHROMOSOMES if chrom not in ["X", "Y", "MT"]],
            ext=wildcards.ext),
        chromosomes = config["resources"] + "ref_fna/chromosomes.list",
    output:
        bcf = config["results"] + "{path}.autosomal.{ext}",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools concat {input.bcfs} \
            -o {output.bcf} \
            -Ob \
        """

# rule vcfgzs_to_autosomal_vcfgz:
#     """Concatenate autosomes."""
#     wildcard_constraints:
#         ext = "vcf.gz",
#     input:
#         vcfs = lambda wildcards: expand(config["results"] + "{path}.chr{chr}.{ext}",
#             path=wildcards.path,
#             chr=[chrom for chrom in CHROMOSOMES if chrom not in ["X", "Y", "MT"]],
#             ext=wildcards.ext),
#         chromosomes = config["resources"] + "ref_fna/chromosomes.list",
#     output:
#         vcf = config["results"] + "{path}.autosomal.{ext}",
#     threads: 1
#     resources: nodes = 1
#     conda: "../envs/common.yaml"
#     shell: """
#         bcftools concat {input.vcfs} \
#             -o {output.vcf} \
#             -Oz \
#         """

rule bcfs_to_autosomal_vcfgz:
    """Concatenate autosomes."""
    input:
        bcfs = lambda wildcards: expand(config["results"] + "{path}.chr{chr}.bcf",
            path=wildcards.path,
            chr=[chrom for chrom in CHROMOSOMES if chrom not in ["X", "Y", "MT"]],
            #ext=wildcards.ext
        ),
        chromosomes = config["resources"] + "ref_fna/chromosomes.list",
    output:
        vcf = config["results"] + "{path}.autosomal.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools concat {input.bcfs} \
            -o {output.vcf} \
            -Oz \
        """

#ruleorder: bcfs_to_autosomal_vcfgz > vcfgzs_to_autosomal_vcfgz

## Ambiguous with biallelics_by_mode
# rule split_chromosomes:
#     """Split into chromosomes."""
#     input:
#         vcf = config["results"] + "{path}/{dataset}.{mode}.vcf.gz",
#         tbi = config["results"] + "{path}/{dataset}.{mode}.vcf.gz.tbi",
#     output:
#         vcf = config["results"] + "{path}/{dataset}.{mode}.chr{chr}.vcf.gz",
#     threads: 1
#     resources: nodes = 1
#     conda: "../envs/common.yaml"
#     shell: """
#         bcftools view {input.vcf} \
#             -r {wildcards.chr} \
#             -Oz \
#             -o {output.vcf}
#         """

# rule number_chromosomes_of_assembly:
#     """Convert chromosomes names to numbered (e.g., 1, 2, 3...) and X, Y, and MT."""
#     input:
#         ref_fasta = config["ref_fasta"],
#         map = config["resources"] + "ref_fna/chr_map.tsv",
#     output:
#         config["resources"] + "ref_fna/numbered.fa.gz",
#     run:
#         # Will result in the command
#         # "gunzip -c [ref_fasta] | sed -e s/[old]/[new]/ -e ... | bgzip > [output]"
#         command = f"gunzip -c {input.ref_fasta} | sed "
#         with open(input.map, 'r') as f:
#             for line in f.readlines():
#                 old, new = line.strip().split("\t")
#                 command += f" -e 's/^>{old}/>{new}/'"
#         command += f" | bgzip > {output}"
#         shell(command)

# rule add_labels_to_pedigree:
#     """Replace animal ids with those including sequence type. For examples, 12345 -> GBS12345 or 11111 -> WGS11111.
#     This doesn't work for when an animal id appears multiple times from different sequencing methods."""
#     input:
#         vcf = config["results"] + "genotypes/subsets/GBS_WES_WGS_labels.284_and_5_swapped.SNP.chr11.min_AF0.01.min_AC2.vcf.gz",
#         pedigree = config["resources"] + "pedigree/{pedigree}.txt",
#     output:
#         pedigree = config["resources"] + "pedigree/{pedigree}.labels.txt",
#     run:
#         samples = shell(f"bcftools query -l {input.vcf}", iterable=True)
#         with open(input.pedigree) as f:
#             text = f.read()
#         for sample in samples:
#             text = text.replace(sample[3:], sample)
#         with open(output.pedigree, "x") as f:
#             f.write(text)




# ## Genetic map






# rule groupby_parents:
#     """Groupby parents with a third column containing offspring separated by commas."""
#     input:
#         #pedigree = config["pedigree"], #config["resources"] + "pedigree/pedigree.2022-12-02_16-13-09.no_underscore.tsv",
#         pedigree = config["results"] + "haplotypes/pedigree/all_samples.all_with_seq.tsv",
#         samples_WGS = config["resources"] + "samples/"
#     output:
#         grouped = config["resources"] + "pedigree/grouped_pedigree.tsv",
#     run:
#         import polars as pl

#         # with open("WGS_no_seq_type.list", "r") as f:
#         #     indivs = f.read().split("\n")

#         data = pl.scan_csv(input.pedigree, separator="\t", has_header=False, infer_schema_length=10000, new_columns=["Indiv", "Sire", "Dam"])

#         data.groupby("Sire", "Dam").agg("Indiv").select(
#             "Sire",
#             "Dam",
#             pl.col("Indiv").arr.join(","),
#         ).collect().write_csv(output.grouped, separator="\t")

# rule family_vcf:
#     input:
#         grouped = config["resources"] + "pedigree/grouped_pedigree.tsv",
#     output:
#         family_vcf = config["results"] + "genetic_map/"
#     shell: """
#         sort -k 2,3 {input.tsv} | 
#         """

#         """
#         read -r SIRE DAM OFFSPRING <<< $()
#         bcftools view {input.vcf} \
#             -s $SIRE,$DAM,$OFFSPRING \
#         """



# rule make_gnotate:
#     """Make gnotate file for SLIVAR."""
#     input:
#         vcfs = "",
#     shell: """
#         slivar make-gnotate \
#             --field AF_popmax:gnomad_popmax_af \
#             --field nhomalt:gnomad_num_homalt \
#             {input.vcfs} \
#         """

# rule slivar:
#     """Use SLIVAR."""
#     input:
#         gnomad = "gnomad_af.zip",
#         ped = config["ped"],  # Pedigree with six columns: fam, sample, pat, mat, sex, and phenotype
#         vcf = "",
#     output:
#         vcf = "",
#     shell: """
#         slivar expr \
#             --gnotate {input.gnomad} \
#             --out-vcf {output.vcf} \
#             --ped {input.ped} \
#         """


rule test:
    """A rule to test if Snakemake is working."""
    output:
        #bcf = config["results"] + "genotypes/pass/{dataset}.12_Indian_12_Chinese_common_to_{subset}.{mode}.chr{chr}.bcf",
        # bcf = config["results"] + "genotypes/pass/{dataset}.{subset}_12_Indian_12_Chinese.{mode}.chr{chr}.bcf",
        test = config["results"] + "test/test.txt",
    threads: 2
    resources:
        nodes = 2
    conda: "../envs/common.yaml"
    shell: """
        echo "This is a test to find the path that `bcftools` is coming from." > {output.test};
        echo $(which bcftools) > {output.test}
        """

rule download_demographics:
    """Pull data from labkey.

    Note that this requires a .netrc file with proper credentials
    in the user home directory in order to connect to the database.
    """
    # params:
    #     species = "Marmoset"  # Will have to rework for baboons that consist of multiple species
    output:
        #config["resources"] + "pedigree/{species}.tsv",
        config["resources"] + "pedigree/demographics.tsv",
    run:
        # This script targets the Python client API version 2.0.0 and later
        import labkey
        from labkey.api_wrapper import APIWrapper
        import pandas as pd
        import polars as pl
        import polars.selectors as cs

        api = APIWrapper("vger.txbiomed.org", "SNPRC", use_ssl=True)

        my_results = api.query.select_rows(
            schema_name="study",
            query_name="Demographics",
            view_name="All Animals",
            columns="Id,species,species/scientific_name,gender,calculated_status,dam,sire,birth,death",
            # filter_array=[
            #     labkey.query.QueryFilter("species/arc_species_code/common_name/value", "Marmoset", "in")
            # ],
            sort="Id"
        )

        df = pl.from_pandas(pd.DataFrame(my_results['rows'])).select(
            ~cs.starts_with("_")
        ).rename({
            # Default names aren't good or may not be consistently names and so are changed here
            "species/scientific_name": "Species",
            "gender": "Sex",
            "birth": "Date of Birth",
            "death": "Date of Death",
            "dam": "Dam",
            "sire": "Sire",
            "calculated_status": "Status",
        }).write_csv(wildcards.output, separator='\t')