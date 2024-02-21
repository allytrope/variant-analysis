rule create_chromosomal_intervals:
    """Create intervals for `rvtest`."""
    input:
        gtf3 = config["gtf3"]
    output:
        set_file = config["resources"] + "annotations/gtf_sets.tsv",
    shell: """
        zcat {input.gtf3} \
            | awk '$3=="gene"' \
            | cut -f 9 \
            | cut -d ';' -f 1 \
            | cut -d ' ' -f 2 \
            | sed 's/"//g' \
            > gtf_ids.list;
        zcat {input.gtf3} \
            | awk '$3=="gene"' \
            | awk 'BEGIN {{OFS=""}} {{print $1,":",$4,"-",$5}}' \
            > gtf_pos.list;
        paste gtf_ids.list gtf_pos.list \
            > {output.set_file};
        """

# For some reason, this rule's `sed` substitutions aren't working
rule create_partial_pheno:
    """Begin creation of phenotype file for `rvtest`."""
    input:
        pops = config["pops"],  # A list of files
        demographics = config["demographics"],
    output:
        pheno = temp(config["resources"] + "samples/pops/partial_pheno.tsv"),
    shell: """
        cat {input.demographics} \
            | cut -f 1,4,6,7 \
            | awk 'BEGIN {{FS="\t"; OFS="\t"}} {{print "F1",$1,$4,$3,$2}}' \
            | sed 's/Male/1/;s/Female/2/;s/Unknown/0/;' \
            | sed '1d' \
            | sed 's/\t\t/\t0\t/g;s/\t\t/\t0\t/g'\
            > {output.pheno};
        """
        # while read LINE; do
        #     indiv=$(echo $LINE | cut -f 1)
        #     if 

        # for FILE in partial.tsv; do \
        # done;

rule create_complete_pheno:
    """Set phenotype by population. Finish phenotype file for `rvtest`."""
    input:
        partial = config["resources"] + "samples/pops/partial_pheno.tsv",
        pops = config["pops"],
    output:
        complete = config["resources"] + "samples/pops/complete_pheno.tsv",
    run:
        with open(input.partial, 'r') as f, open(input.pops[0], 'r') as p1, open(input.pops[1], 'r') as p2, open(output.complete, 'w') as out:
            pop1 = [line.strip() for line in p1.readlines()]
            pop2 = [line.strip() for line in p1.readlines()]
            out.write(" ".join(["fid", "iid", "fatid", "matid", "sex", "y1"]) + "\n")
            for line in f.readlines():
                line_list = line.strip().split("\t")
                indiv = line_list[1]
                if "WGS" + indiv in pop1:
                    pheno = 1
                elif "WES" + indiv in pop1:
                    pheno = 1
                elif "WGS" + indiv in pop2:
                    pheno = 2
                elif "WES" + indiv in pop2:
                    pheno = 2
                else:
                    pheno = 0
                # for seq_type in ["WGS", "WES"]:
                #     line_seq = line_list[:]
                #     line_seq[1] = seq_type + line_list[1]
                #     out.write(" ".join(line_seq + [str(pheno)]) + "\n")
                out.write(" ".join(line_list + [str(pheno)]) + "\n")

rule largest_samples_without_seq_prefix:
    """Create list of largest samples, but remove prefix (by cutting out first three characters)."""
    input:
        samples = config["results"] + "haplotypes/pedigree/{dataset}.largest_samples.list",
    output:
        samples_no_prefix = config["results"] + "haplotypes/pedigree/{dataset}.largest_samples.no_prefix.list",
    threads: 1
    resources: nodes = 1
    conda: "../../envs/bio.yaml"
    shell: """
        cut {input.samples} -c 4- > {output.samples_no_prefix}; \
        """
        # """
        # TMPFILE=$(mktemp); \
        # cut {input.samples} -c 4- > $TMPFILE; \
        # paste {input.samples} $TMPFILE -d ' '> {output.samples_map}; \
        # """
        # """
        # paste {input.samples} $(cut {input.samples} -c 4-) > {output.samples_map}
        # """

rule largest_samples_only_vcf:
    """Remove smaller seq samples. And then rename to remove seq label."""
    input:
        samples = config["results"] + "haplotypes/pedigree/{dataset}.largest_samples.list",
        samples_no_prefix = config["results"] + "haplotypes/pedigree/{dataset}.largest_samples.no_prefix.list",
        vcf = config["results"] + "haplotypes/SHAPEIT5_{seq}/{dataset}.SNP.autosomal.vcf.gz",
        idx = config["results"] + "haplotypes/SHAPEIT5_{seq}/{dataset}.SNP.autosomal.vcf.gz.tbi",
    output:
        vcf = config["results"] + "haplotypes/SHAPEIT5_{seq}/{dataset}.SNP.largest_no_prefix.autosomal.vcf.gz",
    threads: 1
    resources: nodes = 1
    conda: "../../envs/bio.yaml"
    shell: """
        bcftools view {input.vcf} \
            -S {input.samples} \
            -Oz \
        | bcftools reheader \
            -s {input.samples_no_prefix} \
            -o {output.vcf} \
        """


rule kinship_matrix:
    """Create kinship matrix for use in `rvtest`."""
    input:
        bcf = config["results"] + "haplotypes/SHAPEIT5_{seq}/{dataset}.SNP.largest_no_prefix.autosomal.vcf.gz",
        idx = config["results"] + "haplotypes/SHAPEIT5_{seq}/{dataset}.SNP.largest_no_prefix.autosomal.vcf.gz.tbi",
        #ped = config["resources"] + "pheno/RPL.pheno",
        ped = config["ped"],
    params:
        prefix = config["results"] + "stats/burden_test/kinship/{dataset}.{seq}",
    output:
        config["results"] + "stats/burden_test/kinship/{dataset}.{seq}.kinship",
    conda: "../../envs/rvtests.yaml"
    threads: 8
    resources: nodes = 8
    # --bn \
    shell: """
        vcf2kinship \
            --inVcf {input.bcf} \
            --pedigree {input.ped} \
            --thread {threads} \
            --out {params.prefix} \
        """
        #--inVcf {input.vcf} \
        #--bn \

rule burden_test:
    """Note that `test` can be: cmc, zeggini, mb, fp, exactCMC, cmcWald, rarecover, cmat, famcmc, or famzeggini."""
    input:
        vcf = config["results"] + "haplotypes/SHAPEIT5_{seq}/{dataset}.SNP.largest_no_prefix.autosomal.vcf.gz",
        idx = config["results"] + "haplotypes/SHAPEIT5_{seq}/{dataset}.SNP.largest_no_prefix.autosomal.vcf.gz.tbi",
        kinship = config["results"] + "stats/burden_test/kinship/{dataset}.{seq}.kinship",
        pheno = config["resources"] + "pheno/RPL.gt2_pregnancies.pheno",
        set_file = config["resources"] + "annotations/gtf_sets.tsv",
    params:
        out = config["results"] + "stats/burden_test/SHAPEIT5_{seq}/{dataset}.{mode}.burden_test.{test}",
    output:
        config["results"] + "stats/burden_test/SHAPEIT5_{seq}/{dataset}.{mode}.burden_test.{test}.assoc",
        config["results"] + "stats/burden_test/SHAPEIT5_{seq}/{dataset}.{mode}.burden_test.{test}.log",
    conda: "../../envs/rvtests.yaml"
    shell: """
        rvtest \
            --burden {wildcards.test} \
            --inVcf {input.vcf} \
            --kinship {input.kinship} \
            --out {params.out} \
            --pheno {input.pheno} \
            --pheno-name RPL_bool \
            --setFile {input.set_file} \
        """

rule skatO_model:
    """Find SKAT-O. Kernel can be `skat`, `skato`, `kbac`, or `famSkat`."""
    input:
        vcf = config["results"] + "haplotypes/SHAPEIT5_{seq}/{dataset}.SNP.largest_no_prefix.autosomal.vcf.gz",
        idx = config["results"] + "haplotypes/SHAPEIT5_{seq}/{dataset}.SNP.largest_no_prefix.autosomal.vcf.gz.tbi",
        kinship = config["results"] + "stats/burden_test/kinship/{dataset}.{seq}.kinship",
        #pheno = config["resources"] + "samples/pops/complete_pheno.tsv",
        pheno = config["resources"] + "pheno/RPL.gt2_pregnancies.pheno",
        set_file = config["resources"] + "annotations/gtf_sets.tsv",
    params:
        out = config["results"] + "stats/kernel_model/SHAPEIT5_{seq}/{dataset}.{mode}.kernel_model.{test}",
    output:
        config["results"] + "stats/kernel_model/SHAPEIT5_{seq}/{dataset}.{mode}.kernel_model.{test}.assoc",
        config["results"] + "stats/kernel_model/SHAPEIT5_{seq}/{dataset}.{mode}.kernel_model.{test}.log",
    conda: "../../envs/rvtests.yaml"
    shell: """
        rvtest \
            --inVcf {input.vcf} \
            --kernel {wildcards.test} \
            --kinship {input.kinship} \
            --out {params.out} \
            --pheno {input.pheno} \
            --pheno-name RPL_bool \
            --setFile {input.set_file} \
        """

rule single_variant_model:
    """Do single-variant test. Kernel can be `famLRT` or others."""
    input:
        #covar = config["resources"] + "covar/RPL_females.gt2_total_pregnancies.covar",
        vcf = config["results"] + "haplotypes/SHAPEIT5_{seq}/{dataset}.SNP.largest_no_prefix.autosomal.vcf.gz",
        idx = config["results"] + "haplotypes/SHAPEIT5_{seq}/{dataset}.SNP.largest_no_prefix.autosomal.vcf.gz.tbi",
        kinship = config["results"] + "stats/burden_test/kinship/{dataset}.{seq}.kinship",
        pheno = config["resources"] + "pheno/RPL.gt2_pregnancies.pheno",
        set_file = config["resources"] + "annotations/gtf_sets.tsv",
    params:
        out = config["results"] + "stats/single_variant_test/SHAPEIT5_{seq}/{dataset}.{mode}.single_variant.{test}",
    output:
        config["results"] + "stats/single_variant_test/SHAPEIT5_{seq}/{dataset}.{mode}.single_variant.{test}.assoc",
        config["results"] + "stats/single_variant_test/SHAPEIT5_{seq}/{dataset}.{mode}.single_variant.{test}.log",
    conda: "../../envs/rvtests.yaml"
    #            --covar {input.covar} \
    #        --covar-name total_pregnancies,age_at_1st_birth,age_at_last_birth,average_IBI \
    shell: """
        rvtest \
            --covar {input.pheno} \
            --covar-name total_pregnancies \
            --inVcf {input.vcf} \
            --kinship {input.kinship} \
            --out {params.out} \
            --pheno {input.pheno} \
            --pheno-name mortality_rate \
            --setFile {input.set_file} \
            --single {wildcards.test} \
        """


# rule full_ped:
#     """Make full PED from demographics file."""
#     shell: """
#         csvtk -t cut -f 1,7,6,4 {input.demographics} \
#         """

#     run:
#         df = pl.read_csv('Demographics_2023-11-09_16-05-35.tsv', separator='\t', columns=['Id', 'Sire', 'Dam', 'Sex'], dtypes={'Id': str, 'Dam': str, 'Sire': str})        
#         df = df.rename({'Id': 'iid', 'Sex': 'sex', 'Dam': 'matid', 'Sire': 'fatid'})
#         df = df.with_columns(pl.lit('F1').alias('fid')).select(['fid', 'iid', 'fatid', 'matid', 'sex'])
#         df.with_columns(
#             pl.col('sex').str.replace('Female', '2').str.replace('Male', '1').str.replace('Unknown', '0'),
#             pl.col('fatid').fill_null(0),
#             pl.col('matid').fill_null(0),
#         )
#         df.write_csv('from_polars.ped', separator='\t')