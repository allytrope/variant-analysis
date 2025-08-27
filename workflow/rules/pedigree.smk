"""
Rules related to pedigrees
"""


rule SHAPEIT5_fam_file:
    """Create .fam file for SHAPEIT5. This assumes sample names are exactly the animal IDs. And missing fields are filled with 'NA'."""
    input: 
        demographics = config["resources"] + "pedigree/demographics.tsv",
    output:
        fam = config["results"] + "pedigree/all_individuals.fam"
    shell: """
        csvtk cut {input.demographics} \
            -t \
            -f Id,Sire,Dam \
        | csvtk replace \
            -t \
            -f Sire,Dam \
            -p '^$' \
            -r 'NA' \
        > {output.fam} \
        """

rule plink_ped:
    """Create PLINK .ped file from demographics file. All families set to `F1` and all phenotypes set to `0`, which means missing."""
    input:
        demographics = config["resources"] + "pedigree/demographics.tsv",
    output:
        ped = config["resources"] + "pedigree/plink.ped",
    conda: "../envs/common.yaml"
    shell: """
        csvtk cut {input.demographics} \
            -f Id,Sire,Dam,Sex \
            -t \
        | csvtk replace \
            -f Sire,Dam \
            -p '^$' \
            -r '0' \
            -t \
        | sed '1d' \
        | awk 'BEGIN {{FS="\t"; OFS="\t"}} {{print "F1",$0,"0"}}' \
        > {output.ped} \
        """

### Pedigree reconstruction using PRIMUS

rule PLINK_IBD:
    """Create IBD file using PLINK for PRIMUS."""
    input:
        pruned_bcf = config["results"] + "genotypes/pruned/{dataset}.{subset}.{mode}.autosomal.bcf",
    params:
        out_prefix = config["resources"] + "ibd/{dataset}.{subset}.{mode}",
    output:
        ibd = config["results"] + "ibd/{dataset}.{subset}.{mode}.genome",
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
            pl.col("Date of Birth").str.to_date("%m-%d-%Y"),
            today = pl.lit(today).str.to_date('%Y-%m-%d')
        ).with_columns(
            AGE = (pl.col("today") - pl.col("Date of Birth")).dt.total_days(),
            FID = pl.lit("1"),
        ).rename({"Id": "IID"}).select(
            "FID", "IID", "AGE"
        ).write_csv(output.ages, separator="\t")

rule primus:
    """Run PRIMUS."""
    input:
        ibd = config["results"] + "ibd/{dataset}.{subset}.{mode}.genome",
        ages = config["results"] + "pedigree/PRIMUS.ages.tsv",
        sexes = config["results"] + "pedigree/PRIMUS.sexes.tsv",
    output:
        ersa_model = config["results"] + "pedigree_reconstruction/{dataset}.{subset}.{mode}.ERSA_model.txt",
        ersa_results = config["results"] + "pedigree_reconstruction/{dataset}.{subset}.{mode}.ERSA_results.txt",
    #conda: "../envs/primus.yaml"  # TODO: Make this environment
    shell: """
        run_PRIMUS.pl \
            --plink_ibd {input.ibd} \
            --age_file {input.ages} \
            --sex_file {input.sexes} \
            --ersa_model_output {output.ersa_model} \
            --ersa_results {output.ersa_results} \
            --no_IMUS \
        """


# rule padre:
#     """Run PADRE."""


### Pedigree reconstruction with FRANz
# Unfinished
rule create_FRANz_input_file:
    """Get genotypes and birth and death years in the input format for FRANz."""
    input: pruned_bcf = config["results"] + "genotypes/pruned/{dataset}.{subset}.{mode}.autosomal.bcf",
    shell: """
        
        NUM_OF_LOCI

        """

# Unfinished
rule FRANz_pedigree_reconstruction:
    """Reconstruct pedigree using FRANz."""
    input:
        ped = "", # TODO: Add
    params:
        fem_reproductive_age = "6:25", # TODO: Change
        male_reproductive_age = "6:25", # TODO: Change
    shell: """
        FRANz \
            --femrepro {params.fem_reproductive_age} \
            --malerepro {params.male_reproductive_age} \
            --pedigreein {input.ped} \
            --seed 89 \
        """