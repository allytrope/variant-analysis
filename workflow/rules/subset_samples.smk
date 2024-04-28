"""Rules for subsetting samples from VCFs."""

rule list_breeders_and_potential_breeders:
    """Create TSV listing all animals with offspring as well as those that are reserved to be breeders."""
    input:
        demographics = config["demographics"],
    output:
        tsv = config["resources"] + "pedigree/Demographics_2024-04-17_10-28-20.breeders.tsv",
    run:
        import polars as pl

        data = pl.read_csv(input.demographics, separator="\t", infer_schema_length=None).with_columns(
            pl.col("Account Description").str.contains("eserved for breeding").alias("is_reserved_for_breeding"),
            pl.col("Account Description").str.contains("reeder").alias("is_breeder"),
        )
        reserved_for_breeding = data.filter((pl.col("is_reserved_for_breeding") == True)).filter(pl.col("Date of Death").is_null()).select("Id").with_columns(pl.lit(0).cast(pl.UInt32).alias("Num of offspring"))
        sires = data.join(data, left_on="Id", right_on="Sire").group_by("Id").agg(pl.len().alias("Num of offspring"))
        dam = data.join(data, left_on="Id", right_on="Dam").group_by("Id").agg(pl.len().alias("Num of offspring"))
        breeders = pl.concat([sires, dam]).sort("Id")
        pl.concat([breeders, reserved_for_breeding]).sort("Id").write_csv(output.tsv)

rule subset_breeders:
    """Subset out only animals in breeding."""
    input:
        breeders = config["resources"] + "pedigree/Demographics_2024-04-17_10-28-20.breeders.tsv",
        bcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.bcf",
        tbi = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.bcf.csi",
    output:
        bcf = config["results"] + "genotypes/pass/{dataset}_breeders.{mode}.chr{chr}.bcf",
    shell: """
        bcftools view {input.bcf} \
            --force-samples \
            -s $(echo $(cut {input.breeders} -f 1 | sed '1d') | sed 's/ /,/g') \
            -Ob \
            -o {output.bcf} \
        """

rule passing_WGS:
    """Subset out WGS samples."""
    input:
        bcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.bcf",
        csi = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.bcf.csi",
    output:
        bcf = temp(config["results"] + "genotypes/pass/WGS/{dataset}.{mode}.chr{chr}.bcf"),
    shell: """
        bcftools view {input.bcf} \
            -S <(bcftools query -l {input.bcf} | grep WGS) \
            -Ob \
            -o {output.bcf} \
        """

rule passing_WES:
    """Subset out WES samples and keep WES regions determined by mosdepth."""
    input:
        bcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.bcf",
        csi = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.bcf.csi",
        bed = config["results"] + "coverage/common_WES_0.5_loci.bed",  # TODO: Generalize this
    output:
        bcf = temp(config["results"] + "genotypes/pass/WES/{dataset}.{mode}.chr{chr}.bcf"),
    shell: """
        bcftools view {input.bcf} \
            -S <(bcftools query -l {input.bcf} | grep WES) \
            -R {input.bed} \
            -Ob \
            -o {output.bcf} \
        """

rule largest_seq_per_organism:
    """Create a list of samples where only one sample from each organism is kept.
    Only the sequencing type with the most data is kept, which is determined by WGS > WES > GBS > AMP.
    This ordering happens to work here because of the alphabetical ordering.
    Used for determining which sample to use as parent in later step."""
    output:
        largest_samples = config["results"] + "haplotypes/pedigree/{dataset}.largest_samples.list",
    # sed separates id and sequence type. E.g. WGS12345 -> WGS    12345
    # awk flips columns
    # sort rows
    # merge rows based on first column
    # take largest id with largest sequence type (this works because "GBS", "WES", "WGS" are in alphabetical order)
    # grep to remove sample names with an underscore
    # sort
    params:
        samples = '\n'.join(SAMPLES),
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    # sed 's/./&\t/3' {input.samples} \
    shell: """
        echo {params.samples} \
        | sed 's/./&\t/3' \
        | awk -v OFS='\t' '{{print $2,$1}}' \
        | sort \
        | awk '$1!=p{{if(p)print s; p=$1; s=$0; next}}{{sub(p,x); s=s $0}} END {{print s}}' \
        | awk -v OFS='' '{{print $NF,$1}}' \
        | grep -v "_" \
        | sort > {output.largest_samples} \
        """

rule add_seq_to_children:
    """Add seq type (AMP, GBS, WES, and/or WGS) to individual ids in pedigree
    and repeat entries if there are multiple sequencing types for an individual."""
    input:
        #parents = config["pedigree"],
        parents = config["demographics"],
        largest_samples = config["results"] + "haplotypes/pedigree/{dataset}.largest_samples.list",
    output:
        parents = temp(config["results"] + "haplotypes/pedigree/{dataset}.children_with_seq.tsv"),
    params:
        samples = '\n'.join(SAMPLES),
    threads: 1
    resources: nodes = 1
    shell: """
        echo {params.samples} \
        | sed 's/AMP/AMP\t/g;s/GBS/GBS\t/g; s/WES/WES\t/g; s/WGS/WGS\t/g' \
        | sort -k 2 \
        | join - {input.parents} -1 2 -2 1 \
        | awk '{{print $2$1"\t"$3"\t"$4}}' \
        > {output.parents} \
        """
    
rule add_seq_to_parents:
    """Prepend sequence type to parents in pedigree.
    Uses WGS if available, otherwise tries the same sequence type as child.
    If the parent is not sequenced, no sequence type is prepended."""
    input:
        parents = config["results"] + "haplotypes/pedigree/{dataset}.children_with_seq.tsv",
    output:
        tsv = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.tsv",
    params:
        samples = SAMPLES,
    threads: 1
    resources: nodes = 1
    run:
        with open(input.parents, "r") as f:
            with open(output.tsv, "a") as out:
                for line in f:
                    child, sire, dam = line.strip("\n").split("\t")
                    seq_type = child[:3]
                    if f"WGS{sire}" in samples:
                        sire = f"WGS{sire}"
                    elif f"{seq_type}{sire}" in samples:
                        sire = f"{seq_type}{sire}"
                    else:
                        sire = ""
                    if f"WGS{dam}" in samples:
                        dam = f"WGS{dam}"
                    elif f"{seq_type}{dam}" in samples:
                        dam = f"{seq_type}{dam}"
                    else:
                        dam = ""
                    out.write(f"{child}\t{sire}\t{dam}\n")

rule make_forced_ped_format:
    """Add fields to make the correct number of columns for a PLINK PED file.
    However, these extra fields don't actually hold any relevant information.
    They are just a requirement for WhatsHap."""
    input:
        tsv = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.tsv",
    output:
        ped = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.ped",
    threads: 1
    resources: nodes = 1
    shell: """
        awk 'BEGIN {{OFS="\t"}} {{print 0,$1,$2,$3,0,0}}' {input.tsv} \
        | sed 's/\t\t/\t0\t/g' \
        | sed 's/\t\t/\t0\t/g' \
        > {output.ped} \
        """

rule trios_only_tsv:
    input:
        tsv = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.tsv",
    output:
        tsv = config["results"] + "haplotypes/pedigree/{dataset}.trios_only.tsv",
    threads: 1
    resources: nodes = 1
    shell: """
        grep -E "(WGS.*|WES.*|GBS.*|AMP.*){{3}}" {input.tsv} \
        > {output.tsv} \
        """

rule make_trio_vcf:
    """Keep only one child-sire-dam trio."""
    input:
        vcf = config["results"] + "genotypes/pass/{dataset}.{mode}.chr{chr}.vcf.gz",
        tsv = config["results"] + "haplotypes/pedigree/{dataset}.all_with_seq.tsv",
    output:
        trio = temp(config["results"] + "genotypes/pass/trios/chr{chr}/{dataset}.{sample}.{mode}.chr{chr}.vcf.gz"),
    threads: 1
    resources: nodes = 1
    conda: "../envs/bio.yaml"
    shell: """
        read -r CHILD SIRE DAM <<< $(grep ^{wildcards.sample} {input.tsv}); \
        if [ -n "$SIRE" ]; \
        then SIRE=",$SIRE"; \
        fi; \
        if [ -n "$DAM" ]; \
        then DAM=",$DAM"; \
        fi; \
        bcftools view {input.vcf} \
            -s {wildcards.sample}$SIRE$DAM \
            -Oz \
            -o {output.trio} \
        """