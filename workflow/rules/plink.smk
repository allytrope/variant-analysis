"""Rules for PLINK file creation and manipulation. That is, .ped/.map and .bed/.bim/.fam files."""

rule plink_update_parents:
    """For updating parent information in PLINK files."""
    input:
        demographics = config["demographics"],
    output:
        update_parents = temp("{path}/plink/{dataset}.{subset}.update_parents.tsv"),
    conda: "../envs/common.yaml"
    shell: """
        cat {input.demographics} \
        | csvtk cut -f Id,Sire,Dam -t \
        | awk 'BEGIN {{FS="\t"; OFS="\t"}} NR!=1 {{print "F1",$0}} NR==1 {{print "Family",$0}}' \
        | sed 's/\t\t/\t0\t/g' \
        | sed 's/\t$/\t0/g' \
        > {output.update_parents}; \
        """

rule plink_update_sex:
    """For updating sex information in PLINK files."""
    input:
        demographics = config["demographics"],
    output:
        update_sex = temp("{path}/plink/{dataset}.{subset}.update_sex.tsv"),
    conda: "../envs/common.yaml"
    shell:"""
        cat {input.demographics} \
        | csvtk cut -f Id,Sex -t \
        | awk 'BEGIN {{FS="\t"; OFS="\t"}} NR!=1 {{print "F1",$0}} NR==1 {{print "Family",$0}}' \
        | sed 's/Male/1/g;s/Female/2/g;s/Unknown/0/g' \
        > {output.update_sex}; \
        """

rule create_plink_files:
    """Create PLINK .ped and .map files."""
    # `contig` can be chr# or "autosomal"
    input:
        bcf = "{path}/{dataset}.{subset}.{mode}.{contig}.bcf",
        update_parents = "{path}/plink/{dataset}.{subset}.update_parents.tsv",
        update_sex = "{path}/plink/{dataset}.{subset}.update_sex.tsv",
    output:
        ped = "{path}/plink/{dataset}.{subset}.{mode}.{contig}.ped",
        map = "{path}/plink/{dataset}.{subset}.{mode}.{contig}.map",
    params:
        out_prefix = subpath(output.ped, strip_suffix=".ped"),
        ped_or_bed = "",  # Empty string will default to .ped/.map output
    conda: "../envs/rvtests.yaml"
    # --chr {wildcards.chr} \
    shell: """
        plink \
            --bcf {input.bcf} \
            --const-fid F1 \
            --recode \
            {params.ped_or_bed} \
            --out {params.out_prefix} \
            --update-parents {input.update_parents} \
            --update-sex {input.update_sex} \
            --allow-extra-chr \
        """
use rule create_plink_files as create_binary_plink_files with:
    output:
        bed = "{path}/plink/{dataset}.{subset}.{mode}.{contig}.bed",
        bim = "{path}/plink/{dataset}.{subset}.{mode}.{contig}.bim",
        fam = "{path}/plink/{dataset}.{subset}.{mode}.{contig}.fam",
    params:
        out_prefix = subpath(output.bed, strip_suffix=".bed"),
        ped_or_bed = "--make-bed",

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
        | csvtk replace \s
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

# No longer needed since I'm not using an IMPUTE file
# rule impute_to_ped:
#     """Use the `impute_to_ped` command that comes bundled with GERMLINE.
#     This is import to do because using PLINK will arbitrarily change the order of haplotypes
#     according to the GERMLINE documentation."""
#     input:
#         vcf = config["results"] + "haplotypes/Beagle4/{dataset}.{subset}.{mode}.chr{chr}.vcf.gz",
#         samples = config["results"] + "haplotypes/Beagle4/{dataset}.{subset}.{mode}.chr{chr}.samples",
#     params:
#         out_prefix = config["results"] + "haplotypes/Beagle4/{dataset}.{subset}.{mode}.chr{chr}",
#     output:
#         ped = config["results"] + "haplotypes/Beagle4/{dataset}.{subset}.{mode}.chr{chr}.ped",
#         map = config["results"] + "haplotypes/Beagle4/{dataset}.{subset}.{mode}.chr{chr}.map",
#     shell: """
#         impute_to_ped {input.vcf} {input.samples} {params.out_prefix} \
#     """