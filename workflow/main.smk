"""Helpful functions, variables, wildcard constraints, and `include`s for rules."""

import sys
sys.path.insert(0, "./config")
from collections import defaultdict
import os
from pathlib import Path

import polars as pl

# Create defaultdict of runs for samples.
# RUNS = defaultdict(set)
# for file in os.listdir(config["reads"]):
#     split = file.split(".")
#     RUNS[split[0]].add(split[1])


# Read reads table
df = pl.read_csv(config["runs"], comment_prefix='#', separator='\t', schema_overrides={'indiv': pl.String})

def tokenize_sample(fmt):
    """Pull out all values in braces.
    For example, the string `{batch}/{seq}{indiv}_{tissue}` becomes `['batch', 'seq', 'indiv', 'tissue']`.
    """

def wildcardize(fmt):
    """Prepend each variable with 'wildcard.'. This is useful for embedding in shell.
    For sample, '{seq}{indiv}_{library}' becomes '{wildcards.seq}{wildcards.indiv}_{wildcards.library}'.
    """

def collect_samples(fmt, col=None, val=None):
    """Collect into a list values such as sample names and run ids.
    Relies on global `df` variable, which is a polars table
    
    Parameters
    ----------
    fmt : string
        Format in which to output. For example,
        using all fields from table: "{batch}/{seq}{indiv}_{library}_{flowcell_lane}"
    col : string, optional
        Name of a column in the table to filter on
    val : string, optional
        String to filter on
    """
    
    # Parse fmt. This parse words inside braces (including the braces themselves) and strings between reversed braces (the that between "}" and "{"
    tokens = re.findall(r"\{\w+\}|(?<=})[^({|})]+(?={)", fmt)

    polarized_tokens = []
    # Mark tokens as variables or literals
    for token in tokens:
        if "{" in token:
            token = token.replace("{", "").replace("}", "")
            polarized_tokens.append(pl.col(token))
        else:
            polarized_tokens.append(pl.lit(token))
    
    if col == None:
        concatenated_df = df.with_columns(
            concatenated = pl.concat_str(polarized_tokens)
        )
    else:
        # Filter and modify dataframe (uses global df)
        concatenated_df = df.filter(
            pl.col(col) == val
        ).with_columns(
            concatenated = pl.concat_str(polarized_tokens)
        )

    # Keep only unique values and return list
    return list(concatenated_df["concatenated"].unique())


# #Create list of samples
# if config["group_sample_runs_by_batch"]:
#     tmp = 2
# else:
#     tmp = 1
# SAMPLES = sorted(list(set([('_').join(run.split('/')[-1].split('_')[0:tmp]) for run in SAMPLE_RUNS])))

# TODO
# intervals = pl.read_csv(config["resources"] + "ref_fna/contig_intervals.tsv",
#     separator= "\t", columns=["contig", "start", "end"],
#     schema_overrides={"contig": pl.String, "start": pl.Int64, "end": pl.Int64}).with_columns(
#         interval = pl.concat_str([pl.col("contig"), pl.lit(":"), pl.col("start"), pl.lit("-"), pl.col("end")]
#     ))

demographics = pl.read_csv(config["demographics"],
    separator= "\t", columns=["Id", "Sex"], schema_overrides={"Id": pl.String})

# For commands that use bash specific syntax, these have to be activated. Although they may only work with the bio environment.
#shell.executable("/bin/bash")
#shell.prefix("source ~/.bash_profile; ")

wildcard_constraints:
    chr = r"[0-9XYMT]+",
    dataset = r"[A-Za-z0-9\+_\.-]+", #r"\w+",
    filter_method = "hard_filtered|VQSR",
    indiv = r"[A-Za-z0-9\-]+",
    interval = r"[]",
    mode = "SNP|indel|both",
    name = r"[A-Za-z0-9_\.\-\+]+",
    prefix = r"[A-Za-z0-9_/\-\+]+",
    path = r"[A-Za-z0-9\+_/\.-]+",
    #run = r"[A-Za-z0-9_-]+",
    #sample = r"\w+",
    sample1 = r"\w+",
    sample2 = r"\w+",
    seq = "WGS|lpWGS|WES|GBS|AMP|merged",
    start = r"[0-9]+",
    end = r"[0-9]+",
    #library = r"[A-Za-z0-9]+_[A-Za-z0-9]+",
    library = r"[A-Za-z0-9-]+",
    workspace=r"[A-Za-z0-9_\.]+",

#include: "rules/directory_structure.smk"
include: "shortcuts.smk"
include: "rules/sra.smk"
#include: "rules/functions.smk"
include: "rules/indices.smk"
include: "rules/quality_control.smk"
include: "rules/align.smk"
include: "rules/variant_calling.smk"
#include: "rules/variant_recalibration.smk"
include: "rules/hard_filter.smk"
include: "rules/pedigree_reconstruction.smk"
#include: "rules/sra.smk"
#include: "rules/octopus.smk"
include: "rules/phasing.smk"
include: "rules/compression/compress_fastqs.smk"
include: "rules/compression/compress_vcfs.smk"
include: "rules/samples/subset_samples.smk"
include: "rules/samples/superset_samples.smk"
include: "rules/plink.smk"
include: "rules/imputation.smk"
include: "rules/annotation.smk"
include: "rules/liftover.smk"
#include: "rules/gwas.smk"
#include: "rules/gene_counts.smk"
include: "rules/inbreeding.smk"
include: "rules/kinship.smk"
#include: "rules/laser.smk"
include: "rules/relations.smk"
include: "rules/stats/TajimaD.smk"
include: "rules/stats/rare_variant_tests.smk"
include: "rules/scikit-allel.smk"
include: "rules/fst.smk"
include: "rules/pca.smk"
include: "rules/homozygosity.smk"
include: "rules/admixture.smk"
#include: "rules/structural_variation.smk"
#include: "rules/plot_alignments.smk"
#include: "rules/concordance.smk"
#include: "rules/coverage.smk"
#include: "rules/one_off.smk"
include: "rules/misc.smk"

#ruleorder: concat_autosomes > bcf_to_vcfgz
#ruleorder: vep > bcf_to_vcfgz
#ruleorder: vcfgz_to_bcf > call_SVs
#ruleorder: cut_adapters_with_i5_i7_genozipped > cut_adapters_with_i5_i7

# import polars as pl
# srr_ids = list(pl.read_csv("/master/abagwell/variant-analysis/resources/rhesus/samples/X202SC22011650-Z01-F001.srr.list", separator='\t')["srr"])

# rule all:
#     """Generate target files. This rule runs automatically when Snakefile is run.
#     Add or modify the paths under `input` to change the target files generated by snakemake.
#     """
#     input:
#         # Input for generating README.md file
#         config["target_files"],
#         #config["results"] + "README.md",
#         #pass_vcf,
#         #expand(config["resources"] + "subpop/{colony}.samples.list", colony=["rh_SPF_U42", "rh_P51"]),        
#         # expand(config["results"] + "genotypes/pass/{dataset}.{subset}.SNP.dedup_samples.list",
#         #     dataset=config['dataset'],
#         #     subset=config['subset'])
#         # expand(config["results"] + "genotypes/annotated/bcsq/{dataset}.{subset}.SNP.chr{chr}.tsv",
#         #     dataset=config['dataset'],
#         #     subset=config['subset'],
#         #     chr=AUTOSOMES)
#         #king_relatedness_PMx



        
        