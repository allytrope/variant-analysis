import numpy as np
import sgkit as sg
import xarray as xr

ds = sg.load_dataset(snakemake.input.zarr)
cohorts = np.load(snakemake.input.cohorts)
ds["sample_cohort"] = xr.DataArray(cohorts, dims="samples")

ds = sg.window_by_position(ds, size=int(snakemake.wildcards.window_size) * 1000, step=int(snakemake.wildcards.window_size)*500)
#sg.divergence(ds)[["variant_contig", "variant_position", "stat_divergence"]].to_zarr(snakemake.output.zarr)
sg.divergence(ds)["stat_divergence"].to_dataframe().to_pickle(snakemake.output.pickle)