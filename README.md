

## Dependencies

To use all available tools, the following needs to be installed:
* `bwa`
* `conda`
* `mamba`
* `gatk`
* `snakemake`

The method used for testing this program involved installing `conda` and creating a conda environment containing `bwa`, `mamba`, and `snakemake`. To receive the most up-to-date version of `gatk`, download the release from the [GitHub GATK page](https://github.com/broadinstitute/gatk/releases) and build the GATK conda environment as described in [GATK conda instructions](https://gatk.broadinstitute.org/hc/en-us/articles/360035889851--How-to-Install-and-use-Conda-for-GATK4). The environment should be named `envs/gatk.yaml` in the top-level directory of this package.


## Data Preparation

The path and filename of the reference genome should be listed in the appropriate section within `config.yaml`. Similarly, include the path to the input FASTQ files.

Input file paths can be configured in `config.yaml`.

For alignment:
* Reference genome
* Raw FASTQ files (ending with `.fastq.gz`)

For BQSR:
* VCF of known variants for reference genome

For VQSR:
* VCF for truth variants
* VCF for training variantsls


## Output

The output files to be generated (as well as all intermediate files) must be set up under `rule all` in the `Snakefile`. For example, if wanting to create `.vcf` files for each dataset referenced in `config.yaml`, change the string in `rule all` to be `"{output_path}vcf/{dataset}.vcf"`. This requires knowing the directory name as well as the file extension for the desired file. These names can be found under each rule of the `Snakefile` under `output`.


## Running

For running on a SGE cluster, the following command can be issued from the directory that contains the `Snakefile`:
`snakemake --use-conda --cluster "qsub -V -b n -cwd -N smk_variant -o smk_variant.output.log -e smk_variant.error.log" -j 500`

