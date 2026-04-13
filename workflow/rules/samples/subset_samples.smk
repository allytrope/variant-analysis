"""Rules for subsetting samples from VCFs."""

# TODO: Work on this to not conflict with rule samples_names_from_ids
# rule subset_by_colony:
#     """Subset samples by colony."""
#     input:
#         demographics = config["demographics"],
#         runs = config["runs"],
#     output:
#         # "{colony}" can be a value like "rh_SPF_U42" or "rh_P51"
#         samples = config["resources"] + "subpop/samples/{colony}.list",
#     run:
#         import polars as pl

#         # Read colony info
#         demographics = pl.read_csv(input.demographics, separator="\t", columns=["Id", "Colony"], schema_overrides={"Id": pl.String}
#         ).filter(pl.col("Colony") == wildcards.colony)

#         # Read to find unique libraries for individuals
#         libraries = pl.read_csv(input.runs, separator="\t", columns=["indiv", "seq", "library"], schema_overrides={"indiv": pl.String}).unique()

#         # Inner join on the two dataframes
#         demographics.join(libraries, left_on="Id", right_on="indiv", how="inner"
#         ).with_columns(
#             # Construct sample names used in VCF
#             pl.concat_str([pl.col("seq"), pl.col("Id"), pl.lit("_"), pl.col("library")], separator="").alias("sample")
#         ).select("sample"
#         # Write output
#         ).write_csv(output.samples, include_header=False)

rule only_one_seq_type:
    """Subset down to only one sample for each individual and of only a specific sequencing type."""
    wildcard_constraints:
        subset = "AMP|GBS|WES|WGS",
    input:
        vcf = config["results"] + "genotypes/pass/{dataset}.all.{mode}.chr{chr}.bcf",
        csi = config["results"] + "genotypes/pass/{dataset}.all.{mode}.chr{chr}.bcf.csi",
    output:
        vcf = config["results"] + "genotypes/pass/{dataset}.{seq}.{mode}.chr{chr}.bcf",
        tmp_samples = temp(config["results"] + "genotypes/{dataset}.{seq}.{mode}.chr{chr}.samples.list"),
        tmp_samples2 = temp(config["results"] + "genotypes/{dataset}.{seq}.{mode}.chr{chr}.samples2.list"),
    conda: "../../envs/common.yaml"
    shell: """
        bcftools query -l {input.vcf} \
            | sed 's/.\{{3\}}/&\t/; s/_/\t/' \
            | csvtk grep -p {wildcards.seq} -t -f 1 \
            | awk 'BEGIN {{FS="\t"; OFS=""}} $2!=prev {{print $1,$2,"_"$3}} {{prev=$2}}' \
            | tee {output.tmp_samples} \
            | cut -d '_' -f 1 \
            | cut -c 4- \
            > {output.tmp_samples2}; \
        bcftools view {input.vcf} \
            -S {output.tmp_samples} \
        | bcftools reheader  \
            -s {output.tmp_samples2} \
        | bcftools view \
            --min-ac 1 \
            -Ob \
        > {output.vcf} \
        """

rule deduplicate_samples:
    """Create a list of samples where only one samples from each individual is kept."""
    input:
        bcf = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr1.bcf",
        tbi = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr1.bcf.csi",
    output:
        dedup_samples = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.dedup_samples.list",
        dedup_indivs = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.dedup_indivs.list",
    shell: """
            bcftools query -l {input.bcf} \
            | sed 's/.\{{3\}}/&\t/; s/_/\t/' \
            | csvtk sort \
                -k 2,1:u,3:r \
                -L 1:<(for i in WGS lpWGS WES GBS AMP; do echo $i; done) \
                -t \
                -H \
            | awk 'BEGIN {{FS="\t"; OFS=""}} $2!=prev {{print $1,$2,"_"$3}} {{prev=$2}}' \
            | tee {output.dedup_samples} \
            | cut -d '_' -f 1 \
            | cut -c 4- \
            > {output.dedup_indivs}; \
        """

rule deduplicate_individuals:
    """Keep one of each animal and then simplify names down to ids.
    Samples are prioritized in the order: WGS, WES, GBS, AMP.
    Then by the last run name alphabetically (which is usually the newest)."""
    input:
        #bcf = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.SNP.chr{chr}.bcf",
        bcf = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.bcf",
        tbi = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.bcf.csi",
        # vcf = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.vcf.gz",
        # tbi = config["results"] + "genotypes/filtered/{dataset}.{mode}.chr{chr}.vcf.gz.tbi",
    output:
        bcf = config["results"] + "genotypes/pass/{dataset}.{subset}.dedup.{mode}.chr{chr}.bcf",
        tmp_samples = temp(config["results"] + "genotypes/pass/{dataset}.{subset}.dedup.{mode}.chr{chr}.samples.list"),
        tmp_samples2 = temp(config["results"] + "genotypes/pass/{dataset}.{subset}.dedup.{mode}.chr{chr}.samples2.list"),
    conda: "../../envs/common.yaml"
    # Split full sample names into seq + animal_id + run_id assuming a field like so: {seq}{animal_id}_{run_id}
    # Sort by priority
    # Only keep first of each animal_id
    # Then apply list to VCF
    # Reheader sample names to individual names
    shell: """
        bcftools query -l {input.bcf} \
            | sed 's/.\{{3\}}/&\t/; s/_/\t/' \
            | csvtk sort \
                -k 2,1:u,3:r \
                -L 1:<(for i in WGS lpWGS WES GBS AMP; do echo $i; done) \
                -t \
                -H \
            | awk 'BEGIN {{FS="\t"; OFS=""}} $2!=prev {{print $1,$2,"_"$3}} {{prev=$2}}' \
            | tee {output.tmp_samples} \
            | cut -d '_' -f 1 \
            | cut -c 4- \
            > {output.tmp_samples2}; \
        bcftools view {input.bcf} \
            -S {output.tmp_samples} \
        | bcftools reheader  \
            -s {output.tmp_samples2} \
        | bcftools view \
            -Ob \
        > {output.bcf} \
        """
rule match_indivs_to_library:
    """Match a list of indivs to list of libraries.
    This makes it easier to subset VCFs when given a list of individuals
    while the VCF using samples names like {seq}{indiv}_{library}."""
    input:
        indivs = config["resources"] + "subpop/{subset}.indivs.list",
        runs = config["runs"],
    output:
        libraries = config["resources"] + "subpop/{subset}.libraries.list",
    run:
        import polars as pl

        indivs = pl.read_csv(input.indivs, has_header=False, new_columns=["indiv"], separator="\t", comment_prefix="#")
        runs = pl.read_csv(input.runs, separator="\t", comment_prefix="#")
        indivs.join(runs, on="indiv").group_by(
            "indiv"
        ).agg(pl.all().first()).with_columns(
            sample = pl.concat_str([
                pl.col("seq"),
                pl.col("indiv"),
                pl.lit("_"),
                pl.col("library")
            ])
        ).select("sample"
        ).write_csv(output.libraries ,separator="\t", include_header=False)





# rule list_breeders_and_potential_breeders:
#     """Create TSV listing all animals with offspring as well as those that are reserved to be breeders."""
#     input:
#         demographics = config["resources"] + "pedigree/demographics.tsv",
#     output:
#         tsv = config["resources"] + "pedigree/Demographics_2024-04-17_10-28-20.breeders.tsv",
#     run:
#         import polars as pl

#         data = pl.read_csv(input.demographics, separator="\t", infer_schema_length=None).with_columns(
#             pl.col("Account Description").str.contains("eserved for breeding").alias("is_reserved_for_breeding"),
#             pl.col("Account Description").str.contains("reeder").alias("is_breeder"),
#         )
#         reserved_for_breeding = data.filter((pl.col("is_reserved_for_breeding") == True)).filter(pl.col("Date of Death").is_null()).select("Id").with_columns(pl.lit(0).cast(pl.UInt32).alias("Num of offspring"))
#         sires = data.join(data, left_on="Id", right_on="Sire").group_by("Id").agg(pl.len().alias("Num of offspring"))
#         dam = data.join(data, left_on="Id", right_on="Dam").group_by("Id").agg(pl.len().alias("Num of offspring"))
#         breeders = pl.concat([sires, dam]).sort("Id")
#         pl.concat([breeders, reserved_for_breeding]).sort("Id").write_csv(output.tsv)

# rule subset_breeders:
#     """Subset out only animals in breeding."""
#     input:
#         breeders = config["resources"] + "pedigree/Demographics_2024-04-17_10-28-20.breeders.tsv",
#         bcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.bcf",
#         tbi = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.bcf.csi",
#     output:
#         bcf = config["results"] + "genotypes/pass/{dataset}_breeders.{mode}.chr{chr}.bcf",
#     shell: """
#         bcftools view {input.bcf} \
#             --force-samples \
#             -s $(echo $(cut {input.breeders} -f 1 | sed '1d') | sed 's/ /,/g') \
#             -Ob \
#             -o {output.bcf} \
#         """

# rule passing_WES:
#     """Subset out WES samples and keep WES regions determined by mosdepth."""
#     input:
#         bcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.bcf",
#         csi = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.bcf.csi",
#         bed = config["results"] + "coverage/common_WES_0.5_loci.bed",  # TODO: Generalize this
#     output:
#         bcf = temp(config["results"] + "genotypes/pass/WES/{dataset}.{mode}.chr{chr}.bcf"),
#     shell: """
#         bcftools view {input.bcf} \
#             -S <(bcftools query -l {input.bcf} | grep WES) \
#             -R {input.bed} \
#             -Ob \
#             -o {output.bcf} \
#         """

# rule largest_seq_per_organism:
#     """Create a list of samples where only one sample from each organism is kept.
#     Only the sequencing type with the most data is kept, which is determined by WGS > WES > GBS > AMP.
#     This ordering happens to work here because of the alphabetical ordering.
#     Used for determining which sample to use as parent in later step."""
#     output:
#         largest_samples = config["results"] + "haplotypes/pedigree/{dataset}.largest_samples.list",
#     # sed separates id and sequence type. E.g. WGS12345 -> WGS    12345
#     # awk flips columns
#     # sort rows
#     # merge rows based on first column
#     # take largest id with largest sequence type (this works because "GBS", "WES", "WGS" are in alphabetical order)
#     # grep to remove sample names with an underscore
#     # sort
#     params:
#         samples = '\n'.join(collect_samples("{indiv}")),
#     threads: 1
#     resources: nodes = 1
#     conda: "../envs/common.yaml"
#     # sed 's/./&\t/3' {input.samples} \
#     shell: """
#         echo {params.samples} \
#         | sed 's/./&\t/3' \
#         | awk -v OFS='\t' '{{print $2,$1}}' \
#         | sort \
#         | awk '$1!=p{{if(p)print s; p=$1; s=$0; next}}{{sub(p,x); s=s $0}} END {{print s}}' \
#         | awk -v OFS='' '{{print $NF,$1}}' \
#         | grep -v "_" \
#         | sort > {output.largest_samples} \
#         """

# rule add_seq_to_children:
#     """Add seq type (AMP, GBS, WES, and/or WGS) to individual ids in pedigree
#     and repeat entries if there are multiple sequencing types for an individual."""
#     input:
#         #parents = config["pedigree"],
#         parents = config["resources"] + "pedigree/demographics.tsv",
#         largest_samples = config["results"] + "haplotypes/pedigree/{dataset}.largest_samples.list",
#     output:
#         parents = temp(config["results"] + "haplotypes/pedigree/{dataset}.children_with_seq.tsv"),
#     params:
#         samples = '\n'.join(collect_samples("{seq}{indiv}")),
#     threads: 1
#     resources: nodes = 1
#     shell: """
#         echo {params.samples} \
#         | sed 's/AMP/AMP\t/g;s/GBS/GBS\t/g; s/WES/WES\t/g; s/WGS/WGS\t/g' \
#         | sort -k 2 \
#         | join - {input.parents} -1 2 -2 1 \
#         | awk '{{print $2$1"\t"$3"\t"$4}}' \
#         > {output.parents} \
#         """
    
# rule add_seq_to_parents:
#     """Prepend sequence type to parents in pedigree.
#     Uses WGS if available, otherwise tries the same sequence type as child.
#     If the parent is not sequenced, no sequence type is prepended."""
#     input:
#         parents = config["results"] + "haplotypes/pedigree/{dataset}.children_with_seq.tsv",
#     output:
#         tsv = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.tsv",
#     params:
#         samples = collect_samples("{seq}{indiv}"),
#     threads: 1
#     resources: nodes = 1
#     run:
#         with open(input.parents, "r") as f:
#             with open(output.tsv, "a") as out:
#                 for line in f:
#                     child, sire, dam = line.strip("\n").split("\t")
#                     seq_type = child[:3]
#                     if f"WGS{sire}" in samples:
#                         sire = f"WGS{sire}"
#                     elif f"{seq_type}{sire}" in samples:
#                         sire = f"{seq_type}{sire}"
#                     else:
#                         sire = ""
#                     if f"WGS{dam}" in samples:
#                         dam = f"WGS{dam}"
#                     elif f"{seq_type}{dam}" in samples:
#                         dam = f"{seq_type}{dam}"
#                     else:
#                         dam = ""
#                     out.write(f"{child}\t{sire}\t{dam}\n")

# rule make_forced_ped_format:
#     """Add fields to make the correct number of columns for a PLINK PED file.
#     However, these extra fields don't actually hold any relevant information.
#     They are just a requirement for WhatsHap."""
#     input:
#         tsv = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.tsv",
#     output:
#         ped = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.ped",
#     threads: 1
#     resources: nodes = 1
#     shell: """
#         awk 'BEGIN {{OFS="\t"}} {{print 0,$1,$2,$3,0,0}}' {input.tsv} \
#         | sed 's/\t\t/\t0\t/g' \
#         | sed 's/\t\t/\t0\t/g' \
#         > {output.ped} \
#         """

# rule trios_only_tsv:
#     input:
#         tsv = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.tsv",
#     output:
#         tsv = config["results"] + "haplotypes/pedigree/{dataset}.trios_only.tsv",
#     threads: 1
#     resources: nodes = 1
#     shell: """
#         grep -E "(WGS.*|WES.*|GBS.*|AMP.*){{3}}" {input.tsv} \
#         > {output.tsv} \
#         """

# rule make_trio_vcf:
#     """Keep only one child-sire-dam trio."""
#     input:
#         vcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
#         tsv = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.tsv",
#     output:
#         trio = temp(config["results"] + "genotypes/pass/trios/chr{chr}/{dataset}.{sample}.{mode}.chr{chr}.vcf.gz"),
#     threads: 1
#     resources: nodes = 1
#     conda: "../envs/common.yaml"
#     shell: """
#         read -r CHILD SIRE DAM <<< $(grep ^{wildcards.sample} {input.tsv}); \
#         if [ -n "$SIRE" ]; \
#         then SIRE=",$SIRE"; \
#         fi; \
#         if [ -n "$DAM" ]; \
#         then DAM=",$DAM"; \
#         fi; \
#         bcftools view {input.vcf} \
#             -s {wildcards.sample}$SIRE$DAM \
#             -Oz \
#             -o {output.trio} \
#         """

# rule RPL_haplotypes:
#     """Select RPL animal for SHAPEIT5 WGS."""
#     input:
#         #bcf = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.{mode}.chr{chr}.bcf",
#         bcf = config["results"] + "gwas/SHAPEIT5_WGS/bcfs/{dataset}.SNP.chr{chr}.bcf",
#         samples = config["results"] + "samples/RPL_females.WGS.list",
#     output:
#     threads: 1
#     resources: nodes = 1
#     conda: "../envs/common.yaml"
#     shell: """
#         bcftools view {input.bcf} \
#             -S {input.samples} \
#         """