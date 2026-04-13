rule gcta_inbreeding:
    input:
        bed = lambda wildcards: expand(config["results"] + "genotypes/{dir}/plink/{dataset}.{subset}.SNP.chr{chr}.{ext}",
            dir=wildcards.dir,
            chr=AUTOSOMES,
            dataset=wildcards.dataset,
            subset=wildcards.subset,
            ext=["bed", "bim", "fam"]
        ),
        autosomes = config["resources"] + "ref_fna/autosomes.list",
    output:
        mbfile = config["results"] + "inbreeding/GCTA/{dir}/{dataset}.{subset}.files.list",
        ibc = config["results"] + "inbreeding/GCTA/{dir}/{dataset}.{subset}.ibc",
    params:
        prefix = expand(config["results"] + "genotypes/{{dir}}/plink/{{dataset}}.{{subset}}.SNP.chr{chr}",
            chr=AUTOSOMES,
        ),  # Expects .fam, .bim, and .bed
        output = config["results"] + "inbreeding/GCTA/{dir}/{dataset}.{subset}",
    threads: 1
    resources:
        nodes = 1,
        mem_mb = 500_000,  #16_000 was too small
    conda: "../envs/rvtests.yaml"
    shell: """
        for i in {params.prefix}; do echo $i >> {output.mbfile}; done; \

        gcta64 \
            --mbfile {output.mbfile} \
            --autosome-num $(tail -n 1 {input.autosomes}) \
            --autosome \
            --ibc \
            --out {params.output} \
        """
