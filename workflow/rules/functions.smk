## Functions useful for a varietry of rules. Particularly, "collect" functions.
## These are being phased out to be replaced with the rule `collect_samples` found in `Snakefile`.


def collect_runs_from_library(wildcards, full=False):
    """Find libraries from same sample.
    Requires that `wildcards.library` be present in the rule using this function."""
    library = wildcards.library
    runs = []
    for full_run in SAMPLE_RUNS:
        # TODO: Make sure the correct part of the string is tested against
        if library in full_run:
            # For example, pulling out "HV3NHDSXX_L1" from "WES10000_CKDN100000001-1A_HV3NHDSXX_L1"
            if full == False:
                run = "_".join(full_run.split("/")[-1].split("_")[2:])
            else:
                run = full_run
            runs.append(run)
    return runs


def collect_runs_from_sample(wildcards):
    """Find runs from same sample.
    Can take sample name as wildcards.sample or as wildcards.seq + wildcards.indiv.id."""
    try:
        sample = wildcards.sample
    except AttributeError:
        sample = wildcards.seq + wildcards.indiv_id
    try:
        batch = wildcards.batch
        sample = batch + "/" + sample
    except:
        pass
    sample_runs = []
    for run in SAMPLE_RUNS:
        # if wildcards.sample == run.split("_")[0]:
        if sample in run:
            sample_runs.append(sample)
    return sample_runs

def collect_sample_libraries(individual):
    samples = []
    for sample in SAMPLES:
        if individual in sample.split("_")[0]:
            samples.append(sample)
    return ",".join(samples)

