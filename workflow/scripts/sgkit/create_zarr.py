from sgkit.io.vcf import vcf_to_zarr
vcf_to_zarr(snakemake.input.vcf, snakemake.output.zarr)
