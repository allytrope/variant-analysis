## Introduction

This project is a snakemake workflow for processing `.fastq.gz` files and downstream analyses. The variant calling is primarily based on the GATK best practices for germline short variant discovery. Later analyses include pairwise relationship estimation and admixture.


## Dependencies

`conda` is a package manager that can install many of the tools for this program.
Once installed, use the below conda commands to obtain the required packages:
```
conda install -c conda-forge mamba
conda install -c bioconda snakemake
```
Some will still have to be installed manually outside of conda. Depending on which part of the workflow is being used, the following may need to be installed as well and must be accessible through `$PATH`:
* `gatk`
* `gtool`
* `lcMLkin`
* `makeScaffold`
* `salmon`

These that aren't available in conda and must be installed using their corresponding download instructions. For example,`lcMLkin` can be cloned from the [GitHub page](https://github.com/COMBINE-lab/maximum-likelihood-relatedness-estimation). 

Other necessary libraries have packages and versions stored in `workflow/envs`. These will be installed automatically in their own conda environment when running this snakemake pipeline.


## Configuration

Configuration settings are set with a `.yaml` file inside the directory `config`. Multiple config files can exist for convenient testing with different datasets. The name of the desired config file for a particular run, however, should be set for `configfile` in `workflow/Snakefile`. The config file should be filled out to include any starting datasets such as `.fastq.gz` files for all individuals to be analyzed and a reference genome for the species.

Initial required files for first steps of variant calling are the following.
* Reference genome
* Raw FASTQ files (ending with `R1.fastq.gz` and `R2.fastq.gz`)

Other subsequent steps may require their own additional user-created files as well, which can all be specified here in the config file.


## Output

In addition to setting files and variables inside the config files, target outputs (the files that we want to generate) must be listed in the `Snakefile`. When we run, Snakemake will create the files (as well as any intermediates) based on what is listed there.

### Setting Target Files

How to run requires finding the rule in `workflow/rules` for what files are wanted and copying the `output`. For instance, if I want to make .bam files, I can find the corresponding `rule align`, which is in `variant_calling.smk` and find:

```
output:
    alignment = config["results"] + "alignments/raw/{sample}.bam",
```


Now, to place inside list of target files back inside our `Snakefile`, making sure to replace `{sample}` with the actual sample name that we are interested in.

```
rule all:
    input:
        config["results"] + "alignments/raw/WGS12345.bam",
```

Often, many samples will be worked on at the same time. To facilitate this, one could just add more to the target list like so:

```
rule all:
    input:
        config["results"] + "alignments/raw/WGS12345.bam",
        config["results"] + "alignments/raw/WGS23456.bam",
```

Although a more efficient way is to utilize the `expand` command:

```
rule all:
    input:
        expand(config["results"] + "alignments/raw/{sample}.bam", sample=[12345, 23456]),
```

When using many (or all samples), list them in the corresponding file referenced at `config["samples"]` one per line. The file is required for later steps anyway, so helpful to make:
```
rule all:
    input:
        expand(config["results"] + "alignments/raw/{sample}.bam", sample=SAMPLE_NAMES),
```


up until creation of a GenomicsDB datastore in the directory `db/` and subsequent joint calling of variants. 


## Running

Running `snakemake` directly within the project's main directory on the command line will initiate the program by calling `workflow/Snakefile`. To perform a "dry-run", you can use the command `snakemake -n`. This will check that the rules are functional without actually processing any data. This is helpful to quickly test for any obvious problems with new rules. If an upstream file has been modified after a downstream one, Snakemake tries to redo the pipeline from the modified file. If done accidentally, this will delete downstream files, making the rules between run again. So it is good practice to perform the dry-run before each actual run to make sure important files aren't deleted. If a file was accidentally modified and downstream files don't want to be deleted, use Unix command `touch` on all downstream files in order of their creation. This will update the "last modified" timestamps.


### Cluster

To actually run the program on an SGE cluster, the following command can be issued from the project directory:
```
NAME=smk_variant; LOG=log/dirname; nohup snakemake --use-conda --cluster "qsub -V -b n -cwd -pe smp {threads} -N $NAME -o $NAME.out.log -e $NAME.err.log" -j 20 > $LOG/$NAME.smk.log 2>&1
```

The integer in `-j 20` will determine the number of jobs submitted at once to the cluster. Hence a higher number will usually finish quicker, but could take up nodes for others who might be using the cluster. `threads` is replaced with the the integer from the corresponding rule under `threads:` from the corresponding rule. So this can be left alone within the command itself. Set the names for variables `NAME` and `LOG`. `NAME` will be the name given to the SGE and viewable with the `qstat` command. `NAME` also be used for the log files. `LOG` is then the directory in which the log files will be stored. This can be changed as needed to help organize logs.

Besides the generating the files given in the Snakefile's `rule all`, the above command will also generate log files. Check these when errors occur. There are three logs:
* `.smk.log` tells what jobs have been submitted, what files each is trying to make, and when jobs have completed.
* `.err.log` gives details about the actual jobs themselves. When submitting more than one job at once, this can become cluttered with interleaved messages from different jobs. To better find an error message for a specific job, it might require running only a single job at a time, which can be done by setting `-j 1`.
* `.out.log` prints output messages from the jobs. Once again, these will be interleaved. But in general, there is fewer information in the this file and often not as important as the error log.

Make sure to create the directories in which the log files will be placed prior to running to avoid errors.


### Local

Ideally, use the cluster configuration above. Otherwise, to run locally, use: `nohup snakemake --use-conda -c2`. The integer in `-c2` means that two cores will be used. Most rules list the number of threads used, which can be used here. But unlike the cluster method above, this will need to be manually changed. Higher integers (as long as the system has such) will lead to faster processing. This command will leave the error log inside a file in the current directory labelled `nohup.out`.



### Common Errors

`MissingOutputException` - Usually, can just be rerun and will work. Likely due to system latency. Otherwise, if the same file continues to cause this problem, the rule's output (if any) does not match the file(s) listed under the rule's `output`.


