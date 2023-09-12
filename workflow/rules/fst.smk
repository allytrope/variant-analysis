# Find Fst using scikit-allel
rule scikit_allel_fst:
    """Calculate pairwise Fst between populations within non-overlapping windows.
    Output file stores data from a pandas dataframe to be easily graphed with a Python graphing library like seaborn."""
    input: 
        #vcf = CONFIG["vcf"],
        vcf = config["results"] + "structural_variants/SVs/merged/{dataset}.genotyped.pass.bcf",
    output: #config["results"] + "relatedness/fst/created_scikit_fst.txt",
            #hdf5 = config["results"] + "relatedness/fst/fsts.hdf5",
            #csv = config["results"] + "relatedness/fst/fsts.csv",
            pickle = config["results"] + "relatedness/fst/fsts.pickle",
    conda: "../envs/scikit.yaml"
    #params: subpops = CONFIG["fst"]["pops"],
    script: "../scripts/fst.py"

## --> fst_old.ipynb

# Find Fst using vcftools
rule vcftools_fst:
    """Calcuate Fst values between two populations."""
    input:
        #bcf = config["results"] + "structural_variants/SVs/merged/{dataset}.genotyped.pass.bcf",
        bcf = config["results"] + "haplotypes/SHAPEIT5_merged/merged/{dataset}.{mode}.chr{chr}.bcf",
        pop1 = config["pops"][0],
        pop2 = config["pops"][1],
    params:
        out_prefix = config["results"] + "relatedness/fst/{dataset}.{mode}.chr{chr}",
        size = 100_000,
        step = 50_000,
    output:
        out = config["results"] + "relatedness/fst/{dataset}.{mode}.chr{chr}.weir.fst",
    shell: """
        vcftools \
            --bcf {input.bcf} \
            --weir-fst-pop {input.pop1} \
            --weir-fst-pop {input.pop2} \
            --out {params.out_prefix} \
            --fst-window-size {params.size} \
            --fst-window-step {params.step} \
        """

## --> fst.ipynb