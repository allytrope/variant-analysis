## Genozip compression
# rule genozip_fastq:
#     """Genozip FASTQ. Greatly reduces file size, more so than gzip/bgzip."""
#     input:
#         ref_fasta = config["ref_fasta"],
#         R1 = "{path}/{name}.R1.fastq.gz",
#         R2 = "{path}/{name}.R2.fastq.gz",
#     output:
#         genozip = "{path}/{name}.R1+2.fastq.genozip",
#     threads: 8
#     resources: nodes = 8
#     conda: "../envs/common.yaml"
#     # --best ?
#     # --quiet ?
#     shell: """
#         genozip {input.R1} {input.R2} \
#             --pair \
#             --reference {input.ref_fasta} \
#             --threads {threads} \
#         """

# rule genozip_bam:
#     """Genozip BAM, replacing original. Greatly reduces file size, more so than gzip/bgzip."""
#     input:
#         ref_fasta = config["ref_fasta"],
#         bam = "{path}/{name}.bam",
#     output:
#         bam = "{path}/{name}.bam",
#         bai = "{path}/{name}.bam.bai",
#     threads: 8
#     resources: nodes = 8
#     conda: "../envs/common.yaml"
#     # --best ?
#     # --quiet ?
#     shell: """
#         genozip {input.bam} \
#             --reference {input.ref_fasta} \
#             --replace \
#             --threads {threads} \
#         """

# rule genounzip_bam:
#     """Unzip genozip-compressed BAM."""
#     input:
#         ref_fasta = config["ref_fasta"],
#         geno = "{path}/{name}.bam.geno",
#     output:
#         bam = "{path}/{name}.bam",
#         bai = "{path}/{name}.bam.bai",
#     threads: 8
#     resources: nodes = 8
#     conda: "../envs/common.yaml"
#     shell: """
#         genounzip {input.geno} \
#             --index \
#             --reference {input.ref_fasta} \
#             --threads {threads} \
#         """


rule genounzip_to_fastq_gz:
    """Unzip genozipped R1+2 FASTQ into R1 or R2 as a pipe.
    Note: When using the piped output, sometimes it doesn't work because it takes too long."""
    wildcard_constraints:
        read = "R1|R2",
    input:
        fastq = config["reads"] + "{batch}/{seq}{indiv}_{library}_{flowcell_lane}.R1+2.fastq.genozip",
        ref_genozip = config["compression"]["ref_fasta"],
    output: 
        #fastq = temp(pipe(config["resources"] + "reads/{batch}/{seq}{indiv}_{library}_{flowcell_lane}.{read}.fastq.gz")),
        fastq = temp(config["resources"] + "reads/{batch}/{seq}{indiv}_{library}_{flowcell_lane}.{read}.fastq.gz"),
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        genocat {input.fastq} --{wildcards.read} --reference {input.ref_genozip} | bgzip -c > {output.fastq}
        """
        
rule genounzip_fastq:
    """Unzip genozipped R1+2 FASTQ into R1 or R2 as a pipe."""
    wildcard_constraints:
        read = "R1|R2",
    input:
        fastq = config["reads"] + "{batch}/{seq}{sample_run}.R1+2.fastq.genozip",
        ref_genozip = config["compression"]["ref_fasta"],
    output: 
        fastq = temp(pipe(config["resources"] + "reads/{batch}/{seq}{sample_run}.{read}.fastq")),
        #fastq = config["resources"] + "reads/{batch}/{seq}{sample_run}.{read}.fastq",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        genocat {input.fastq} --{wildcards.read} --reference {input.ref_genozip} > {output.fastq}
        """

rule genounzip_fastq_interleaved:
    """Unzip genozipped R1+2 FASTQ into an interleaved FASTQ as a pipe."""
    wildcard_constraints:
        read = "R1|R2",
    input:
        fastq = config["reads"] + "{batch}/{seq}{sample_run}.R1+2.fastq.genozip",
        ref_genozip = config["compression"]["ref_fasta"],
    output: 
        fastq = pipe(config["resources"] + "reads/{batch}/{seq}{sample_run}.R1+R2.fastq"),
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        genocat {input.fastq} --reference {input.ref_genozip} > {output.fastq}
        """

# rule genounzip_fastq_one_line:
#     """Unzip first line only of genozipped R1+2 FASTQ as a pipe."""
#     wildcard_constraints:
#         read = "R1|R2",
#     input:
#         fastq = config["reads"] + "{batch}/{seq}{sample_run}.R1+2.fastq.genozip",
#         ref_genozip = config["compression"]["ref_fasta"],
#     output: 
#         fastq = temp(config["results"] + "reads/{batch}/{seq}{sample_run}.first_line.R1+2.fastq"),
#     threads: 1
#     resources: nodes = 1
#     conda: "../envs/common.yaml"
#     shell: """
#         genocat {input.fastq} --reference {input.ref_genozip} | head -n 1 > {output.fastq}
#         """

rule fastq_one_line:
    """Unzip first line only of genozipped R1+2 FASTQ as a pipe."""
    wildcard_constraints:
        read = "R1|R2",
    input:
        fastq = config["reads"] + "{batch}/{seq}{sample_run}.R1.fastq.gz",
        ref_genozip = config["compression"]["ref_fasta"],
    output: 
        fastq = temp(config["results"] + "reads/{batch}/{seq}{sample_run}.first_line.R1+2.fastq"),
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        zcat {input.fastq} | head -n 1 > {output.fastq}
        """