# This Snakefile is meant to be run initially when setting up a project.
# Subsequent analyses should be run through workflow/{project}.Snakefile.

rule chromosome_windows:
    """Create sliding windows for chromosomes."""
    input:
        fai = config["ref_fasta"] + ".fai",
    output:
        windows = config["resources"] + "ref_fna/contig_intervals.tsv",
    params:
        window_size = 50_000_000,
    # NOTE: The path doesn't have the "../" prefix when I specify this Snakefile with `-s <Snakefile>`
    conda: "../envs/common.yaml"
    shell: """
        awk -v window={params.window_size} '{{ \
            chr=$1; \
            len=$2; \
            for (start=1; start<len; start+=window) {{ \
                end=start+window-1; \
                if (end>len) end=len; \
                print chr"\\t"start"\\t"end; \
            }} \
        }}' {input.fai} \
        > {output.windows} \
        """

rule all:
    """Initial files to be created before running rest of workflow."""
    input:
        config["resources"] + "ref_fna/contig_intervals.tsv"

    