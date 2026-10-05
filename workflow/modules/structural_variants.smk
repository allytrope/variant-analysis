"""Rules for structural variants calling."""

include: "rules/indices.smk"
# include: "rules/pedigree_reconstruction.smk"

# Structural variants
include: "rules/structural_variants/long-read_SV_calling.smk"
include: "rules/structural_variants/short-read_SV_calling.smk"
include: "rules/structural_variants/merging_SVs.smk"

include: "rules/compression/compress_vcfs.smk"
include: "rules/plink.smk"
include: "rules/phasing.smk"
include: "rules/imputation.smk"
include: "rules/misc.smk"