

## Dependencies

To use all available tools, the following need to be installed and accessible through $PATH:
* snakemake
* fastp
* bwa
* docker

Also install the official GATK image through Docker.
After Docker is installed, this can be done by executing the command `docker pull broadinstitue/gatk`.


## Data Preparation

Reference genome and files to be processed should be listed in their appropriate sections within `workflow/config.yaml`.

Input file paths can be configured in `config.yaml`
* Reference genome -> `data/ref`
* Raw FASTQ files -> `data/fastq` (files should be compressed, ending with `.fq.gz`)


## Running

This data analysis pipeline can be run using the `snakemake` command followed by the file to be generated. For example, `snakemake vcf/specimen.vcf` would generate the VCF file for the dataset "specimen".