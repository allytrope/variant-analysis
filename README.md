

## Dependencies

To use all available tools, the following needs to be installed and accessible through $PATH:
* `bwa`
* `conda`
* `mamba`
* `gatk`
* `snakemake`


## Data Preparation

The path and filename of the reference genome should be listed in the appropriate section within `config.yaml`. Similarly, include the path to the input FASTQ files.

Input file paths can be configured in `config.yaml` for:
* Reference genome
* Raw FASTQ files

Similarly, the output directory can also be specified for the final files and any intermediates.
FASTQ files should also be gzipped and ending with `.fastq.gz`.


## Running

This data analysis pipeline can be run using the `snakemake --use-conda` command followed by the file to be generated. For example, `snakemake ../vcf/specimen.vcf` would generate the VCF file from the dataset "specimen".