"""
Rules related to pedigree reconstruction.
"""


### Pedigree reconstruction using PRIMUS

# Maybe I should just use the indiv IDs
# rule replace_sample_underscores:
#     """Rename samples in VCF by replacing underscores with hyphens.
#     This is necessary for PLINK, which otherwise splits sample names
#     on underscores and uses the first half as FID and the second half as IID."""
#     input:
#         bcf = config["results"] + "genotypes/pruned/{dataset}.{subset}.{mode}.autosomal.bcf",
#     output:
#         bcf = config["results"] + "genotypes/pruned/{dataset}.{subset}.{mode}.autosomal.underscores_removed.bcf",
#     conda: "../envs/common.yaml"
#     shell: """
#         bcftools view {input.bcf} -Ov \
#         | bcftools reheader -s <(bcftools query -l {input.bcf} \
#             | awk 'BEGIN {{OFS="\t"}} {{print $1,$1}}' \
#             | sed 's/_/-/2') \
#         | bcftools view -Ob -o {output.bcf} \
#     """

# Alternative to using prePRIMUS. However, it split sample names on underscores
# with the first half for FID and the second for IID
rule PLINK_IBD:
    """Create IBD file using PLINK for PRIMUS."""
    input:
        pruned_bcf = config["results"] + "genotypes/pruned/{dataset}.{subset}.dedup.{mode}.autosomal.bcf",
    params:
        out_prefix = subpath(output.ibd, strip_suffix=".genome"),
    output:
        ibd = config["results"] + "ibd/{dataset}.{subset}.dedup.{mode}.genome",
    conda: "../envs/common.yaml"
    shell: """
        plink \
            --bcf {input.pruned_bcf} \
            --genome \
            --out {params.out_prefix} \
            --noweb \
            --allow-extra-chr \
        """

rule PRIMUS_sex_file:
    """Generate tables of sexes for each animal for PRIMUS.
    Has columns FID, IID, and SEX. For SEX, 0 = unknown, 1 = male, and 2 = female."""
    input:
        demographics = config["resources"] + "pedigree/demographics.tsv",
    output:
        sexes = config["results"] + "pedigree/PRIMUS.sexes.tsv"
    shell: """
        csvtk cut {input.demographics} \
            -t \
            -f Id,Sex \
        | csvtk replace \
            -t \
            -f Sex \
            -p Unknown \
            -r 0 \
        | csvtk replace \
            -t \
            -f Sex \
            -p Male \
            -r 1 \
        | csvtk replace \
            -t \
            -f Sex \
            -p Female \
            -r 2 \
        | csvtk rename \
            -t \
            -f Id,Sex \
            -n ID,SEX \
        | sed 's/^/1\t/1' \
        | csvtk rename \
            -t \
            -f '1' \
            -n FID \
        > {output} \
    """

rule PRIMUS_age_file:
    """Generate tables of ages for each animal for PRIMUS.
    Has columns FID, IID, and AGE. AGE refers to how old an individual is or would be if it didn't die in days."""
    input:
        demographics = config["resources"] + "pedigree/demographics.tsv",
    output:
        ages = config["results"] + "pedigree/PRIMUS.ages.tsv"
    run:
        from datetime import date
        import polars as pl

        today = str(date.today())
        df = pl.read_csv(input.demographics, separator='\t', columns=["Id", "Date of Birth"], #"Date of Death"],
            schema_overrides={"Id": pl.String}
        ).with_columns(
            pl.col("Date of Birth").str.split(" ").list.get(0).str.to_date("%Y-%m-%d"),  #.to_date("%m-%d-%Y"),
            today = pl.lit(today).str.to_date('%Y-%m-%d')
        ).with_columns(
            AGE = (pl.col("today") - pl.col("Date of Birth")).dt.total_days(),
            FID = pl.lit("1"),
        ).rename({"Id": "IID"}).select(
            "FID", "IID", "AGE"
        ).write_csv(output.ages, separator="\t")

# No longer needed since I'm not using an IMPUTE file
# rule write_samples:
#     """Input sample file for `impute_to_ped` command. Not sure why though it needs this since sample names are in a VCF already
#     Also, `impute_to_ped` skips the first two lines of this file for some reason."""
#     input:
#         vcf = config["results"] + "haplotypes/Beagle5/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz",
#     output:
#         samples = config["results"] + "haplotypes/Beagle5/{dataset}.{subset}.{mode}.chr{chr}.samples",
#     shell: """
#         echo '# `impute_to_ped` removes these first' >> {output.samples};
#         echo '# two lines for unknown reasons' >> {output.samples};
#         bcftools query -l {input.vcf} \
#         | cut -c 4- \
#         | cut -d _ -f 1 \
#         >> {output.samples}
#         """

rule GERMLINE:
    """Run GERMLINE for ERSA."""
    input:
        # TODO: These need to be after SHAPEIT5
        ped = config["results"] + "haplotypes/SHAPEIT5/plink/{dataset}.{subset}.{mode}.chr{chr}.ped",
        map = config["results"] + "haplotypes/SHAPEIT5/plink/{dataset}.{subset}.{mode}.chr{chr}.map",
    output:
        match = config["results"] + "ibd/GERMLINE/{dataset}.{subset}.{mode}.chr{chr}.match",
    params:
        out_prefix = subpath(output.match, strip_suffix=".match"),
        out_dir = config["results"] + "ibd/GERMLINE",
    log:
        config["results"] + "ibd/GERMLINE/{dataset}.{subset}.{mode}.chr{chr}.log",
    # Some values taken from https://github.com/belowlab/tnprc-pedigrees/blob/master/code.sh
    # `germline` won't work unless the output directory already exists, hence the use of `mkdir`
    shell: """
        mkdir -p {params.out_dir}; \
        germline \
            -input {input.ped} {input.map} \
            -min_m 2.5 \
            -err_het 1 \
            -err_hom 2 \
            -output {params.out_prefix} \
    """

rule GERMLINE_segments:
    """Create list of files with GERMLINE segments for all autosomes."""
    input:
        segments = lambda wildcards: expand("{results}ibd/GERMLINE/{dataset}.{subset}.{mode}.chr{chr}.match",
            results=config["results"],
            dataset=wildcards.dataset,
            subset=wildcards.subset,
            mode=wildcards.mode,
            #chr=CHROMOSOMES,
            chr=AUTOSOMES,
            ext=["", ".tbi"],),
    output:
        #temp(config["results"] + "ibd/GERMLINE/{dataset}.{subset}.{mode}.segments.list"),
        config["results"] + "ibd/GERMLINE/{dataset}.{subset}.{mode}.segments.list",
    shell: """
        for FILE in {input.segments}; do \
            echo $FILE >> {output}; \
        done \
        """

# Attempting to use `container: singularity`, but not working
# rule COMPADRE:
#     """Run COMPADRE. This serves to replace PRIMUS, ERSA, and PADRE in one step."""
#     input:
#         ibd = config["results"] + "ibd/{dataset}.{subset}.{mode}.genome",
#         ages = config["results"] + "pedigree/PRIMUS.ages.tsv",
#         sexes = config["results"] + "pedigree/PRIMUS.sexes.tsv",
#         segments = config["results"] + "ibd/GERMLINE/{dataset}.{subset}.{mode}.segments.list",
#     params:
#         output_dir = subpath(output.test, strip_suffix=".test"),
#     output:
#         # TODO: Just a test to see if rule is working
#         test = config["results"] + "pedigree_reconstruction/COMPADRE/{dataset}.{subset}.{mode}.test",
#     #singularity:
#     container:
#         "/master/abagwell/singularity_images/compadre.sif"
#     # Some values taken from https://github.com/belowlab/tnprc-pedigrees/blob/master/code.sh
#     # TODO: Generalize number of chromosomes
#     # --no_IMUS \
#     #         perl run_COMPADRE.pl \
#     shell: """
#             --plink_ibd {input.ibd} \
#             --segment_data {input.segments} \
#             --run_padre \
#             --age_file {input.ages} \
#             --sex_file {input.sexes} \
#             --number_of_chromosomes=20 \
#             --rec_per_meioses=13.6239623 \
#             --confidence_level=0.999 \
#             --output {params.output_dir}; \
#         touch {output.test} \
#         """

# Using singularity inside of `run:`
rule COMPADRE:
    """Run COMPADRE. This serves to replace PRIMUS, ERSA, and PADRE in one step."""
    input:
        ibd = config["results"] + "ibd/{dataset}.{subset}.dedup.{mode}.genome",
        ages = config["results"] + "ibd/PRIMUS.ages.tsv",
        sexes = config["results"] + "ibd/PRIMUS.sexes.tsv",
        segments = config["results"] + "ibd/GERMLINE/{dataset}.{subset}.dedup.{mode}.segments.list",
    params:
        input_dir = '/master/abagwell/variant-analysis/results/rhesus/ibd',
        #output_dir = subpath(output.test, strip_suffix=".test"),
        output_dir = '/master/abagwell/variant-analysis/results/rhesus/pedigree_reconstruction/COMPADRE',
        sif = "/master/abagwell/singularity_images/compadre.sif",
    output:
        # TODO: Just a test to see if rule is working
        test = config["results"] + "pedigree_reconstruction/COMPADRE/{dataset}.{subset}.dedup.{mode}.test",
    #singularity:
    # Some values taken from https://github.com/belowlab/tnprc-pedigrees/blob/master/code.sh
    # TODO: Generalize number of chromosomes
    # --no_IMUS \
    # --genome
    shell: """
        singularity run \
            --bind {params.input_dir}:/compadre_data \
            --bind {params.output_dir}:/compadre_output \
        {params.sif} \
            --segment_data /compadre_data/GERMLINE/{wildcards.dataset}.{wildcards.subset}.dedup.{wildcards.mode}.segments.list \
            --plink_ibd /compadre_data/{wildcards.dataset}.{wildcards.subset}.dedup.{wildcards.mode}.genome \
            --run_padre \
            --age_file /compadre_data/PRIMUS.ages.tsv \
            --sex_file /compadre_data/PRIMUS.sexes.tsv \
            --number_of_chromosomes=20 \
            --rec_per_meioses=13.6239623 \
            --confidence_level=0.999 \
            --output /compadre_output
            --port_number 22; \
        touch {output.test} \
        """
