"""Rules for processing VCFs."""

include: "../rules/hard_filter.smk"
include: "../rules/phasing.smk"
include: "../rules/compression/compress_vcfs.smk"
include: "../rules/plink.smk"