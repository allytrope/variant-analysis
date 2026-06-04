# Find Fst using scikit-allel
rule FSTest_fst:
    """Calculate Fst with FSTest.
    
    Will have to modify the commands depending on what populations are being compared."""
    input:
        bcf = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.bcf",
        demographics = config["resources"] + "pop/MML_groups_from_Martha.fixed4.with_Brooks_origin.tsv",
    params:
    #     pop1 = "bcftools view {input.bcf} -Ov -S <(csvtk grep {input.demographics} -t -f Interval -p Founders | csvtk grep -t -f Origin -p Brooks | csvtk -t cut -f Id | sed '1d');",
    #     pop2 = "bcftools view {input.bcf} -Ov -S <(csvtk grep {input.demographics} -t -f Interval -p Founders | csvtk grep -t -f Origin -p 'non-Brooks' | csvtk -t cut -f Id | sed '1d');",
    #     pop3 = "bcftools view {input.bcf} -Ov -S <(csvtk grep {input.demographics} -t -f Interval -p Founders2 | csvtk -t cut -f Id | sed '1d');",
        out = config["results"] + "fst/FSTest/{dataset}.{subset}.{mode}.chr{chr}.method{method}",
    # Method must be 1 through 4
    output:
        out = config["results"] + "fst/FSTest/{dataset}.{subset}.{mode}.chr{chr}.method{method}.fst",
    conda: "../envs/scikit.yaml"
    shell: """
        python3 ~/bin/FSTest1.3.py \
            --pop1 <(bcftools view {input.bcf} --force-samples -Ov -S <(csvtk grep {input.demographics} -t -f Interval -p Founders | csvtk grep -t -f Origin -p Brooks | csvtk -t cut -f Id | sed '1d')) \
            --pop2 <(bcftools view {input.bcf} --force-samples -Ov -S <(csvtk grep {input.demographics} -t -f Interval -p Founders | csvtk grep -t -f Origin -p 'non-Brooks' | csvtk -t cut -f Id | sed '1d')) \
            --m {wildcards.method} \
            --zt 1 \
            --o {params.out} \
        """


rule scikit_allel_fst:
    """Calculate pairwise Fst between populations within non-overlapping windows.
    Output file stores data from a pandas dataframe to be easily graphed with a Python graphing library like seaborn."""
    input: 
        #vcf = CONFIG["vcf"],
        vcf = config["results"] + "structural_variants/SVs/merged/{dataset}.genotyped.pass.bcf",
    output: #config["results"] + "relatedness/fst/created_scikit_fst.txt",
            #hdf5 = config["results"] + "relatedness/fst/fsts.hdf5",
            #csv = config["results"] + "relatedness/fst/fsts.csv",
            pickle = config["results"] + "fst/allel/fsts.pickle",
    conda: "../envs/scikit.yaml"
    #params: subpops = CONFIG["fst"]["pops"],
    script: "../scripts/fst.py"

## --> fst_old.ipynb

# Find Fst using vcftools
rule vcftools_fst:
    """Calcuate Fst values between two populations.
    For this rule, use BCF that has filtered out on MAF<0.01 or possibly 0.05.
    Although doesn't seem to actually work for more than 2 populations.
    It will make output, but won't show any additional pairwise Fst values."""
    input:
        #bcf = config["results"] + "structural_variants/SVs/merged/{dataset}.genotyped.pass.bcf",
        #bcf = config["results"] + "haplotypes/SHAPEIT5_merged/merged/{dataset}.{mode}.chr{chr}.bcf",
        bcf = config["results"] + "genotypes/pass/{dataset}.{subset}.{mode}.chr{chr}.bcf",
        #pop1 = config["pops"][0],
        #pop2 = config["pops"][1],
        #demographics = config["resources"] + "pop/MML_groups_from_Martha.fixed4.tsv",
        demographics = config["resources"] + "pop/MML_groups_from_Martha.fixed4.with_Brooks_origin.tsv",
    params:
        out_prefix = config["results"] + "fst/{dataset}.{subset}.{mode}.chr{chr}",
        pop1 = config["results"] + "fst/{dataset}.{subset}.{mode}.chr{chr}.pop1.list",
        pop2 = config["results"] + "fst/{dataset}.{subset}.{mode}.chr{chr}.pop2.list",
        pop3 = config["results"] + "fst/{dataset}.{subset}.{mode}.chr{chr}.pop3.list",
        #size = 100_000,
        #step = 50_000,
    output:
        out = config["results"] + "fst/{dataset}.{subset}.{mode}.chr{chr}.weir.fst",
    # The creation of the pop files will change based on what file is used and how it is filtered.
        # csvtk grep {input.demographics} -t -f Interval -p Founders | csvtk grep -t -f Origin -p Brooks | csvtk -t cut -f Id | sed '1d' > {params.pop1}; \
        # csvtk grep {input.demographics} -t -f Interval -p Founders | csvtk grep -t -f Origin -p 'non-Brooks' | csvtk -t cut -f Id | sed '1d' > {params.pop2}; \
        # csvtk grep {input.demographics} -t -f Interval -p Founders2 | csvtk -t cut -f Id | sed '1d' > {params.pop3}; \
    shell: """        
        vcftools \
            --bcf {input.bcf} \
            --weir-fst-pop {params.pop2} \
            --weir-fst-pop {params.pop3} \
            --out {params.out_prefix} \
        """
    #         --fst-window-size {params.size} \
    #         --fst-window-step {params.step} \
    

## --> fst.ipynb
