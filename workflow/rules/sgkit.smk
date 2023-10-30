"""Rules for calculating stats using sgkit."""

rule create_zarr:
    input:
        vcf = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.{mode}.autosomal.bcf",
        csi = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.{mode}.autosomal.bcf.csi",
    output:
        zarr = directory(config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.{mode}.autosomal.zarr"),
    threads: 20
    resources: nodes = 20
    conda: "../envs/sgkit.yaml"
    script:
        "../scripts/sgkit/create_zarr.py"

rule cohorts_npy:
    """Create cohorts file to be used by sgkit. Automatically deletes itself after used."""
    input:
        zarr = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.{mode}.autosomal.zarr",
        pops = config["pops"],
    output:
        #npy = temp(config["resources"] + "samples/pops/{dataset}.{mode}.cohorts.npy"),
        npy = config["resources"] + "samples/pops/{dataset}.{mode}.cohorts.npy",
    conda: "../envs/sgkit.yaml"
    script:
        "../scripts/sgkit/cohorts_npy.py"

# Currently, windowing is not working correctly.
rule sgkit_divergence:
    """Estimate nucleotide divergence between two populations."""
    input:
        zarr = config["results"] + "haplotypes/SHAPEIT5_WGS/{dataset}.{mode}.autosomal.zarr",
        cohorts = config["resources"] + "samples/pops/{dataset}.{mode}.cohorts.npy"
    output:
        pickle = config["results"] + "sgkit/divergence/{dataset}.{mode}.{window_size}kbp.pickle",
    threads: 20
    resources: nodes = 20
    conda: "../envs/sgkit.yaml"
    script:
        "../scripts/sgkit/divergence.py"


rule render_sgkit_divergence:
    input:
        pickle = config["results"] + "sgkit/divergence/{dataset}.{mode}.{window_size}kbp.pickle",
    output:
        svg = config["results"] + "sgkit/divergence/{dataset}.{mode}.{window_size}kbp.svg",
    # conda:
    #     "../envs/graph.yaml"
    notebook:
        "../notebooks/sgkit/divergence.ipynb"