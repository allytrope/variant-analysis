## Dependencies

`conda` is a package manager that can install many of the tools for this program.
Once installed, `conda install`:
* `mamba`
* `snakemake`

Depending on which part of the workflow the workflow is being used, the following may need to be installed as well and must be accessible through `$PATH`:
* `admixture`
* `bcftools`
* `bwa`
* `gatk`
* `gtool`
* `impute2`
* `lcMLkin`
* `plink`
* `shapeit`
* `vcftools`

To receive the most up-to-date version of `gatk`, download the release from the [GitHub GATK page](https://github.com/broadinstitute/gatk/releases). GATK's own dependecies will be installed automatically in their own environment when needed at `/workflow/envs/gatk.yaml`. Some of these other tools can be installed through `conda`. Those that aren't available can be installed using their corresponding download instructions. For example,`lcMLkin` can be cloned from the [GitHub page](https://github.com/COMBINE-lab/maximum-likelihood-relatedness-estimation). 


## Configuration

Configuration settings are set with a `.yaml` file inside the directory `config`. Multiple config files can exist for convenient testing with different datasets. The name of the desired config file for a particular run, however, should be set under `configfile:` in `workflow/Snakefile`. The config file should be filled out to include any starting datasets such as `.fastq` files for all individuals to be analyzed and a reference genome for the species. Other steps may require their own additional files as well and these can all be specified here. 

Some of the required files for various steps include the following.
For alignment:
* Reference genome
* Raw FASTQ files (ending with `.fastq.gz`)

For BQSR:
* VCF of known variants for reference genome

For VQSR:
* VCF for truth variants
* VCF for training variants


## Output

The output files to be generated (as well as all intermediate files) must be set up under `rule all` in the `Snakefile`. For example, if wanting to create `.vcf` files for each dataset referenced in `config.yaml`, change the string in `rule all` to be `"{results}vcf/{dataset}.vcf"`. This requires knowing the directory name as well as the file extension for the desired file. These names can be found under each rule of the `.smk` files in `workflow/rules` under `output`.


## Running

Running `snakemake` directly within the project's main directory on the command line will initiate the program by calling `workflow/Snakefile`. To perform a "dry-run", you can use the command `snakemake -n`. This will check that the rules are functional without actually processing any data. This is can be helpful to quickly test for any obvious problems with new rules.

To actually run the program on a SGE cluster, the following command can be issued from the directory:
`VAR=smk_variant; snakemake --use-conda --cluster "qsub -V -b n -cwd -pe smp {threads} -N $VAR -o $VAR.output.log -e $VAR.error.log" -j 100`

The integer in `-j 100` will determine the number of jobs submitted at once to the cluster. Hence a higher number will usually finish quicker, but could take up nodes for others who might be using the cluster. `threads` is replaced with the the integer from the corresponding rule. `-o` and `-e` refer to the paths to output and error logs respectively. These may also be changed. Check these when errors occur. Often snakemake will just say that a nonzero exit code has been given. In order to find the problem, the snakemake job will have to be run outside of a cluster. The command for executing this way is: `snakemake --use-conda -c2`. The integer in `-c2` refers to the number of threads that this process will use. This can be increased as needed. Since some processes may take a long time, it is suggested that you use `nohup snakemake --use-conda -c2`. With this command, even if a connection is lost, snakemake will continue the command and leave the error log inside a file labelled `nohup.out`.
