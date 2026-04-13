"""Workflow for determing structural variation in genomes, that is, deletions, insertions, inversions, tandem duplications, and translocations."""

# rule bcf_to_vcf:
#     """Stream a VCF from a BCF."""
#     input:
#         bcf = "{path}.bcf"
#     output:
#         vcf = temp(pipe("{path}.vcf"))
#     threads: 1
#     resources:
#         nodes = 1
#     shell: """
#         samtools view {input.bcf} \
#             -Ov \
#             -o {output.vcf} \
#         """



rule truvari:
    """Comparing SVs using Truvari."""
    input:
        base_calls = config["results"] + "structural_variants/delly/merged/LRS.vcf.gz",  # Must end in .gz and be bgzipped
        base_calls_idx = config["results"] + "structural_variants/delly/merged/LRS.vcf.gz.tbi",
        comp_calls = config["results"] + "structural_variants/delly/merged/4WGS.vcf.gz",  # Must end in .gz and be bgzipped
        comp_calls_idx = config["results"] + "structural_variants/delly/merged/4WGS.vcf.gz.tbi",  # Must end in .gz and be bgzipped
    output:
        config["results"] + "structural_variants/test/test.txt",
        # config["results"] + "structural_variants/truvari/tp-base.vcf.gz",
        # config["results"] + "structural_variants/truvari/tp-comp.vcf.gz",
        # config["results"] + "structural_variants/truvari/fp.vcf.gz",
        # config["results"] + "structural_variants/truvari/fn.vcf.gz",
        # config["results"] + "structural_variants/truvari/summary.json",
        #config["results"] + "structural_variants/truvari/params.json",
    params:
        out_dir = config["results"] + "structural_variants/truvari/",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly2.yaml"
    shell: """
        truvari bench \
            -b {input.base_calls} \
            -c {input.comp_calls} \
            -o {params.out_dir} \
        """


# ----------------------#
#  SURVIVOR (consensus) #
# ----------------------#

rule SV_consensus:
    """Find consensus between SV callers."""
    input:
        # list_of_vcfs = lambda wildcards: expand(config["results"] + "structural_variants/{SV_caller}/merged/{dataset}.vcf",
        #     #SV_caller=["delly", "sniffles"],  # For LRS
        #     SV_caller=["delly", "smoove"],  # For WGS (non-LRS)
        #     dataset=wildcards.dataset),
        list_of_vcfs = expand(config["results"] + "structural_variants/delly/merged/{dataset}.vcf", dataset=["LRS", "4WGS"]),
    output:
        merged = config["results"] + "structural_variants/SURVIVOR/{dataset}.vcf",
    params:
        max_distance = 1000,  # in bases
        min_consensus = 2,
        match_type = 0,  # 0=False, 1=True
        match_strand = 1,  # 0=False, 1=True
        unused = 0,  # Current versions of SURVIVOR do not use this parameter anymore
        min_length = 0,  # Minimum length of SVs
    threads: 24
    resources: nodes = 24
    conda: "../envs/delly2.yaml"
    shell: """
        SURVIVOR merge <(for FILE in {input.list_of_vcfs}; do echo $FILE; done) \
            {params.max_distance} \
            {params.min_consensus} \
            {params.match_type} \
            {params.match_strand} \
            {params.unused} \
            {params.min_length} \
            {output.merged} \
        """

# rule all_SURVIVOR_combinations:
#     """Summarize overlap between all callers. Requires running every combination of 2 or more callers with SURVIVOR."""
#     input: "",
#     output: "",
#     threads: 24
#     resources: nodes = 24
#     conda: "../envs/delly2.yaml"
#     shell: """

#         """


# NanoSV
# rule call_SVs_NanoSV:
#     """SV call LRS data using NanoSV."""
#     threads: 8
#     shell: """
#         NanoSV \
#             -t {threads} \
#         """



# ---------#
#  GRIDSS  #
# ---------#

# Under development
rule GRIDSS_setupreference:
    """Call short-read structural variants with GRIDSS"""
    input:
        ref_fasta = "/master/abagwell/variant-analysis/resources/rhesus/ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa"
    output:
        #dict = config["ref_fasta"] + ".dict",
        dict = "/master/abagwell/variant-analysis/resources/rhesus/ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.dict",
    threads: 1
    resources:
        nodes = 1
    conda: "../envs/gridss.yaml"
    shell: """
        gridss \
            --reference {input.ref_fasta} \
            -s setupreference \
        """

rule GRIDSS_preprocess:
    """Call short-read structural variants with GRIDSS"""
    input:
        bams = collect_runs_from_sample,  # From variant_calling.smk
        bams_idx = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_sample(wildcards))),
        #ref_fasta = config["ref_fasta"],  #
        ref_fasta = "/master/abagwell/variant-analysis/resources/rhesus/ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa",  # Must not be gzipped
        dict = "/master/abagwell/variant-analysis/resources/rhesus/ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.dict",
        #dict = config["ref_fasta"] + ".dict",
    output:
        bam = config["results"] + "structural_variants/GRIDSS/{sample}.bam.gridss.working/{sample}.bam.sv.bam",
    params:
        #labels = lambda wildcards, input: ",".join([sample.split("_")[0] for sample in SAMPLES]),
        labels = lambda wildcards, input: ",".join([sample.split("_")[0] for sample in SAMPLE_RUNS]),
        output_dir = config["results"] + "structural_variants/GRIDSS",
    threads: 4
    resources:
        nodes = 4
    conda: "../envs/gridss.yaml"
    # --skipsoftcliprealignment works because the BAMs were aligned with BWA-mem
    shell: """
        gridss {input.bams} \
            --labels {params.labels} \
            --reference {input.ref_fasta} \
            -s preprocess \
            --skipsoftcliprealignment \
            --threads {threads} \
            --workingdir {params.output_dir} \
        """

def split_samples(jobnodes: int):
    """Split samples of size `jobnodes` or `jobnodes + 1` for the remainders."""
    split_list = []
    for i in range(0, len(SAMPLES), jobnodes):  
        split_list.append(SAMPLES[i:i + jobnodes])

    # Copy samples from last node into the previous
    if len(split_list[-1]) < jobnodes:
        for idx, sample in enumerate(split_list[-1]):
            split_list[idx].append(sample)
    
    # Remove last node
    split_list.pop()

    return split_list


def collect_runs_from_sample2(sample):
    """Find runs (files) from same sample."""
    sample_runs = []
    for run in SAMPLE_RUNS:
        # if wildcards.sample == run.split("_")[0]:
        if sample in run:
            sample_runs.append(config["results"] + "alignments/recalibrated/" + run + ".bam")
    return sample_runs

def collect_run_names_from_sample2(sample):
    """Find runs (names) from same sample."""
    sample_runs = []
    for run in SAMPLE_RUNS:
        # if wildcards.sample == run.split("_")[0]:
        if sample in run:
            sample_runs.append(run)
    return sample_runs

def collect_runs_for_group(wildcards):
    runs = []
    for sample in split_samples(int(wildcards.jobnodes))[int(wildcards.jobindex)]:
        #runs.append(*collect_runs_from_sample2(sample))
        runs += collect_runs_from_sample2(sample)
    return runs

def collect_run_names_for_group(wildcards):
    runs = []
    for sample in split_samples(int(wildcards.jobnodes))[int(wildcards.jobindex)]:
        #runs.append(*collect_run_names_from_sample2(sample))
        runs += collect_run_names_from_sample2(sample)
    return runs

rule GRIDSS_assembly:
    """Call short-read structural variants with GRIDSS"""
    input:
        # bams = collect_runs_from_sample,  # From variant_calling.smk
        # bams_idx = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_sample(wildcards))),
        bams = lambda wildcards: collect_runs_for_group(wildcards),
        bam_indices = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_for_group(wildcards))),
        #bams = expand(config["results"] + "alignments/recalibrated/{sample}.bam", sample=SAMPLE_RUNS),
        #bais = expand(config["results"] + "alignments/recalibrated/{sample}.bam.bai", sample=SAMPLE_RUNS),
        preprocess_bams = expand(config["results"] + "structural_variants/GRIDSS/{sample}.bam.gridss.working/{sample}.bam.sv.bam", sample=SAMPLE_RUNS),
        #ref_fasta = config["ref_fasta"],
        ref_fasta = "/master/abagwell/variant-analysis/resources/rhesus/ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa",
        ref_dict = "/master/abagwell/variant-analysis/resources/rhesus/ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.dict",
        ref_fai = "/master/abagwell/variant-analysis/resources/rhesus/ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.fai",
    output:
        bam = config["results"] + "structural_variants/GRIDSS/assembly.{jobindex}of{jobnodes}.bam",
    params:
        #labels = lambda wildcards, input: ",".join([sample.split("_")[0] for sample in SAMPLE_RUNS]),
        labels = lambda wildcards, input: ",".join([sample.split("_")[0] for sample in collect_run_names_for_group(wildcards)]),
        output_dir = config["results"] + "structural_variants/GRIDSS",
    threads: 8
    resources:
        nodes = 8
    conda: "../envs/gridss.yaml"
    # --skipsoftcliprealignment works because the BAMs were aligned with BWA-mem
    shell: """
        gridss {input.bams} \
            --assembly {output.bam} \
            --labels {params.labels} \
            --reference {input.ref_fasta} \
            -s assemble \
            --threads {threads} \
            --workingdir {params.output_dir} \
        """

# rule gather_assemblies:
#     """Gather assembly results together."""
#     input:
#     output:
#     params:
#         output_dir = config["results"] + "structural_variants/GRIDSS",
#     threads: 8
#     resources:
#         nodes = 8
#     conda: "../envs/gridss.yaml"
#     shell: """
#         gridss {input.bams} \
#             -a {input.assembly} \
#             -s assemble \
#             --threads {threads} \
#             --workingdir {params.output_dir} \
#         """

rule merge_BAM_assemblies:
    """Merge batched BAM assemblies for input into GRIDSS."""
    input:
        assembled_bams = lambda wildcards: expand(config["results"] + "structural_variants/GRIDSS/assembly.{jobindex}of{jobnodes}.bam",
            jobindex=range(0, int(wildcards.jobnodes)),
            jobnodes=wildcards.jobnodes),
        assembled_bam_indicies = lambda wildcards: expand(config["results"] + "structural_variants/GRIDSS/assembly.{jobindex}of{jobnodes}.bam.bai",
            jobindex=range(0, int(wildcards.jobnodes)),
            jobnodes=wildcards.jobnodes),
    output:
        merged_bam = config["results"] + "structural_variants/GRIDSS/assembly.merged.{jobnodes}-nodes.bam"
    shell: """
        samtools merge {input.assembled_bams} \
            -o {output.merged_bam} \
        """

rule GRIDSS_calling:
    """Call short-read structural variants with GRIDSS."""
    input:
        # assembled_bams = lambda wildcards: expand(config["results"] + "structural_variants/GRIDSS/assembly.{jobindex}of{jobnodes}.bam",
        #     jobindex=range(0, int(wildcards.jobnodes)),
        #     jobnodes=wildcards.jobnodes),
        # assembled_bam_indicies = lambda wildcards: expand(config["results"] + "structural_variants/GRIDSS/assembly.{jobindex}of{jobnodes}.bam.bai",
        #     jobindex=range(0, int(wildcards.jobnodes)),
        #     jobnodes=wildcards.jobnodes),
        bams = expand(config["results"] + "alignments/recalibrated/{sample}.bam", sample=SAMPLE_RUNS),
        bais = expand(config["results"] + "alignments/recalibrated/{sample}.bam.bai", sample=SAMPLE_RUNS),
        #ref_fasta = config["ref_fasta"],
        merged_bam = config["results"] + "structural_variants/GRIDSS/assembly.merged.{jobnodes}-nodes.bam",
        merged_bam_bai = config["results"] + "structural_variants/GRIDSS/assembly.merged.{jobnodes}-nodes.bam.bai",
        ref_fasta = "/master/abagwell/variant-analysis/resources/rhesus/ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa",
        ref_dic = "/master/abagwell/variant-analysis/resources/rhesus/ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.dict",
    output:
        vcf = config["results"] + "structural_variants/GRIDSS/batched.{jobnodes}-nodes.vcf.gz",
    params:
        #assemblies = lambda wildcards, input: " ".join([f"--assembly {bam.split('/')[-1]}" for bam in input.assembled_bams]),
        labels = lambda wildcards, input: ",".join([sample.split("_")[0] for sample in SAMPLE_RUNS]),
        output_dir = config["results"] + "structural_variants/GRIDSS",
    threads: 8
    resources:
        nodes = 8
    conda: "../envs/gridss.yaml"
    # --skipsoftcliprealignment works because the BAMs were aligned with BWA-mem
    #            {params.assemblies} \
    shell: """
        gridss {input.bams} \
            --assembly {input.merged_bam} \
            --labels {threads} \
            --output {output.vcf} \
            --reference {input.ref_fasta} \
            -s call \
            --threads {threads} \
            --workingdir {params.output_dir} \
        """




# ---------#
#  smoove  #
# ---------#

rule smoove_call:
    """Call structural variants with LUMPY through smoove."""
    input:
        # bams = collect_runs_from_sample,  # From variant_calling.smk
        # bams_idx = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_sample(wildcards))),
        bam = config["results"] + "alignments/merged/{sample}.bam",
        bai = config["results"] + "alignments/merged/{sample}.bam.bai",
        ref_fasta = config["ref_fasta"],
    output:
        #bcf = config["results"] + "structural_variants/smoove/per_sample/{seq}{indiv_id}.bcf",
        vcf = config["results"] + "structural_variants/smoove/per_sample/{sample}-smoove.genotyped.vcf.gz",
        #histo = config["results"] + "structural_variants/smoove/per_sample/{sample}.histo",
    params:
        outdir = config["results"] + "structural_variants/smoove/per_sample/",
    threads: 1
    resources: nodes = 1
    conda: "../envs/smoove_simplified.yaml"
    #container: "docker://brentp/smoove"
    #smoove call {input.bams} \
    # NOTE: Process subsitution does not here. Also had problems adding multiple BAMs separately.
    shell: """
        smoove call {input.bam} \
            --fasta {input.ref_fasta} \
            --genotype \
            --name {wildcards.sample} \
            --outdir {params.outdir} \
            -p {threads} \
        """

rule smoove_merge:
    """Merge samples from `smoove call` runs."""
    input:
        vcfs = expand(config["results"] + "structural_variants/smoove/per_sample/{sample}-smoove.genotyped.vcf.gz", sample=SAMPLES),
        ref_fasta = config["ref_fasta"],
    output:
        #bcf = config["results"] + "structural_variants/smoove/per_sample/{seq}{indiv_id}.bcf",
        vcf = config["results"] + "structural_variants/smoove/merged_ungenotyped/{dataset}.vcf.gz",
    params:
        outdir = config["results"] + "structural_variants/smoove/merged_ungenotyped/",
    threads: 1
    resources: nodes = 1
    conda: "../envs/smoove_simplified.yaml"
    shell: """
        smoove merge {input.vcfs} \
            --fasta {input.ref_fasta} \
            --name merged \
            --outdir {params.outdir}; \
        mv {params.outdir}merged.sites.vcf.gz {params.outdir}{wildcards.dataset}.vcf.gz \
        """

rule smoove_genotype:
    """Merge samples from `smoove call` runs."""
    input:
        # bams = collect_runs_from_sample,  # From variant_calling.smk
        # bams_idx = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_sample(wildcards))),
        bam = config["results"] + "alignments/merged/{sample}.bam",
        bai = config["results"] + "alignments/merged/{sample}.bam.bai",
        ref_fasta = config["ref_fasta"],
        vcf = config["results"] + "structural_variants/smoove/merged_ungenotyped/{dataset}.vcf.gz",
    output:
        #bcf = config["results"] + "structural_variants/smoove/per_sample/{seq}{indiv_id}.bcf",
        vcf = config["results"] + "structural_variants/smoove/all_sites/{dataset}/{sample}-smoove.genotyped.vcf.gz",
    params:
        outdir = config["results"] + "structural_variants/smoove/all_sites/{dataset}",
    threads: 4
    resources: nodes = 4
    conda: "../envs/smoove_simplified.yaml"
    shell: """
        smoove genotype {input.bam} \
            -d \
            --fasta {input.ref_fasta} \
            --name {wildcards.sample} \
            --outdir {params.outdir} \
            -p 4 \
            --vcf {input.vcf} \
            -x \
        """

rule smoove_paste:
    """Merge samples that now have all sites across samples."""
    input:
        vcfs = lambda wildcards: expand(
            config["results"] + "structural_variants/smoove/all_sites/{dataset}/{sample}-smoove.genotyped.vcf.gz",
                dataset=wildcards.dataset,
                sample=SAMPLES),
    output:
        #bcf = config["results"] + "structural_variants/smoove/per_sample/{seq}{indiv_id}.bcf",
        bcf = config["results"] + "structural_variants/smoove/merged/{dataset}.smoove.square.vcf.gz",
    params:
        outdir = config["results"] + "structural_variants/smoove/merged/",
    threads: 4
    resources: nodes = 4
    conda: "../envs/smoove_simplified.yaml"
    shell: """
        smoove paste {input.vcfs} \
            --name {wildcards.dataset} \
            --outdir {params.outdir}; \
        mv {params.outdir}{wildcards.dataset}.smoove.square.vcf.gz {params.outdir}{wildcards.dataset}.vcf.gz \
        """


# Manta
def add_option(bams):
    options = ""
    print("Testing")
    for bam in bams:
        options += f"--bam {bam} "
    return options
rule configure_Manta:
    input:
        bams = expand(config["results"] + "alignments/recalibrated/{sample_run}.bam", sample_run=SAMPLE_RUNS),
        ref_fasta = config["ref_fasta"],
    output:
        command = config["results"] + "structural_variants/manta/runWorkflow.py",
    params:
        bam_options = lambda wildcards, input: add_option(input.bams),
        out_dir = config["results"] + "structural_variants/manta/",
    threads: 12
    resources: nodes = 12
    conda: "../envs/manta.yaml"
    shell: """
        configManta.py \
            {params.bam_options} \
            --referenceFasta {input.ref_fasta} \
            --runDir {params.out_dir} \
        """

rule execute_Manta:
    input:
        command = config["results"] + "structural_variants/manta/runWorkflow.py",
    output:
        diploid_SVs = config["results"] + "structural_variants/manta/results/variants/diploidSV.vcf.gz",
    threads: 8
    resources: nodes = 8
    conda: "../envs/manta.yaml"
    shell: """
        {input.command} \
            -j {threads} \
        """
        


#---------#
#  DELLY  #
#---------#

rule call_SVs:
    """Call SVs per sample (for short-read WGS)."""
    wildcard_constraints:
        seq = "WGS",
    input:
        bams = collect_runs_from_sample,  # From variant_calling.smk
        bams_idx = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_sample(wildcards))),
        ref_fasta = config["ref_fasta"],
    output:
        bcf = config["results"] + "structural_variants/delly/per_sample/{seq}{indiv_id}.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly2.yaml"
    shell: """
        delly call {input.bams} \
            -g {input.ref_fasta} \
            -o {output.bcf} \
        """

rule call_SVs_LRS:
    """Call SVs per sample (for ONT LRS)."""
    wildcard_constraints:
        seq = "LRS",
    input:
        bams = collect_runs_from_sample,  # From variant_calling.smk
        bams_idx = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_sample(wildcards))),
        ref_fasta = config["ref_fasta"],
    output:
        bcf = config["results"] + "structural_variants/delly/per_sample/{seq}{indiv_id}.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly2.yaml"
    shell: """
        delly lr {input.bams} \
            -y ont \
            -g {input.ref_fasta} \
            -o {output.bcf} \
        """

rule merge_SVs:
    """Merge BCFs with structural variants."""
    input:
        bcfs = expand(config["results"] + "structural_variants/delly/per_sample/{sample}.bcf", sample=SAMPLES),
    output:
        bcf = config["results"] + "structural_variants/delly/merged_ungenotyped/{dataset}.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly2.yaml"
    shell: """
        delly merge {input.bcfs} \
            -o {output.bcf} \
        """

rule genotype_merged_SVs:
    """Genotype merged BCFs with structural variants (for short-read WGS)."""
    wildcard_constraints:
        seq = "WGS",
    input:
        bams = collect_runs_from_sample,  # From variant_calling.smk
        bams_idx = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_sample(wildcards))),
        bcf = config["results"] + "structural_variants/delly/merged_ungenotyped/{dataset}.bcf",
        ref_fasta = config["ref_fasta"],
    output:
        bcf = config["results"] + "structural_variants/delly/all_sites/{dataset}.{seq}{indiv_id}.genotyped.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly2.yaml"
    shell: """
        delly call {input.bams} \
            -g {input.ref_fasta} \
            -v {input.bcf} \
            -o {output.bcf} \
        """

rule genotype_merged_SVs_LRS:
    """Genotype merged BCFs with structural variants (for ONT LRS)."""
    wildcard_constraints:
        seq = "LRS",
    input:
        bams = collect_runs_from_sample,  # From variant_calling.smk
        bams_idx = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_sample(wildcards))),
        # bam = config["results"] + "alignments/recalibrated/merged/{seq}{indiv_id}.bam",
        # bai = config["results"] + "alignments/recalibrated/merged/{seq}{indiv_id}.bam.bai",
        bcf = config["results"] + "structural_variants/delly/merged_ungenotyped/{dataset}.bcf",
        ref_fasta = config["ref_fasta"],
    output:
        bcf = config["results"] + "structural_variants/delly/all_sites/{dataset}.{seq}{indiv_id}.genotyped.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly2.yaml"
    shell: """
        delly lr {input.bams} \
            -y ont \
            -g {input.ref_fasta} \
            -v {input.bcf} \
            -o {output.bcf} \
        """

rule merge_all_sites_SVs:
    """Merge samples with sites from all other samples."""
    input:
        bcfs = lambda wildcards: expand(
            config["results"] + "structural_variants/delly/all_sites/{dataset}.{sample}.genotyped.bcf",
                dataset=wildcards.dataset,
                sample=SAMPLES),
    output:
        bcf = config["results"] + "structural_variants/delly/merged/{dataset}.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools merge {input.bcfs}\
            -m id \
            -Ob \
            -o {output.bcf} \
        """

rule filter_SVs:
    """Filter structural variants."""
    input:
        bcf = config["results"] + "structural_variants/delly/merged/{dataset}.bcf",
        csi = config["results"] + "structural_variants/delly/merged/{dataset}.bcf.csi",
    output:
        bcf = config["results"] + "structural_variants/delly/merged/{dataset}.filtered.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly2.yaml"
    shell: """
        delly filter {input.bcf} \
            -f germline \
            -o {output.bcf} \
        """

rule passing_SVs:
    """Filter structural variants."""
    input:
        bcf = config["results"] + "structural_variants/delly/merged/{dataset}.filtered.bcf",
        csi = config["results"] + "structural_variants/delly/merged/{dataset}.filtered.bcf.csi",
    output:
        bcf = config["results"] + "structural_variants/delly/merged/{dataset}.pass.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools view {input.bcf} \
            -f \
            -Ob \
            -o {output.bcf} \
        """

rule split_by_SV_type:
    """Split merged BCF by type of SV."""
    wildcard_constraints:
        SV_type = "DEL|INS|DUP|INV|BND|ALL",
    input:
        bcf = config["results"] + "structural_variants/delly/merged/{dataset}.filtered.bcf",
    output:
        bcf = config["results"] + "structural_variants/delly/merged/{dataset}.{SV_type}.pass.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools view {input.bcf} \
            -i "INFO/SVTYPE='<{wildcards.SV_type}>'" \
            -f \
            -Oz \
            -o {output.bcf} \
        """

## Sniffles

## Using Sniffles
# rule sniffles:
#     """Call SVs for long read sequences on single sample."""
#     wildcard_constraints:
#         SV_type = "DEL|INS|DUP|INV|BND|ALL",
#     input:
#         bam = "/master/abagwell/variant-analysis/resources/rhesus/reads/LRS/{sample}.final.bam",
#         bai = "/master/abagwell/variant-analysis/resources/rhesus/reads/LRS/{sample}.final.bam.bai",
#         ref_fasta = config["ref_fasta"],
#     output:
#         #vcf = config["results"] + "structural_variants/SVs/sniffles/{dataset}.{SV_type}.genotyped.pass.vcf.gz",
#         vcf = config["results"] + "structural_variants/SVs/sniffles/{sample}.vcf.gz",
#     threads: 1
#     resources: nodes = 1
#     conda: "../envs/sniffles.yaml"
#     shell: """
#         sniffles \
#             -i {input.bam} \
#             -v {output.vcf} \
#             --reference {input.ref_fasta} \
        # """

rule sniffles_SNF:
    """Call SVs as SNF file for long read sequences on a sample."""
    input:
        bams = collect_runs_from_sample,  # From variant_calling.smk
        bams_idx = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_sample(wildcards))),
        # bam = "/master/abagwell/variant-analysis/resources/rhesus/reads/LRS/{sample}.final.bam",
        # bai = "/master/abagwell/variant-analysis/resources/rhesus/reads/LRS/{sample}.final.bam.bai",
        ref_fasta = config["ref_fasta"],
    output:
        #vcf = config["results"] + "structural_variants/SVs/sniffles/{dataset}.{SV_type}.genotyped.pass.vcf.gz",
        snf = config["results"] + "structural_variants/sniffles/per_sample/{sample}.snf",
    threads: 24
    resources: nodes = 24
    conda: "../envs/rvtests.yaml"
    shell: """
        sniffles \
            --input <(samtools merge {input.bams} -o -) \
            --snf {output.snf} \
            --reference {input.ref_fasta} \
            --threads {threads} \
        """

rule sniffles_merge_SNFs:
    """Merge SNF files."""
    input:
        snfs = lambda wildcards: expand(
            config["results"] + "structural_variants/sniffles/per_sample/{sample}.snf",
                sample=SAMPLES),
        ref_fasta = config["ref_fasta"],
    output:
        #vcf = config["results"] + "structural_variants/SVs/sniffles/{dataset}.{SV_type}.genotyped.pass.vcf.gz",
        vcf = config["results"] + "structural_variants/sniffles/merged/LRS.bcf",
    threads: 24
    resources: nodes = 24
    conda: "../envs/rvtests.yaml"
    shell: """
        sniffles \
            --input {input.snfs} \
            --vcf >(samtools -Ob -o {output.vcf}) \
            --reference {input.ref_fasta} \
            --threads {threads} \
        """

rule split_sniffles_by_SV_type:
    """Split merged BCF by type of SV."""
    wildcard_constraints:
        SV_type = "DEL|INS|DUP|INV|BND|ALL",
    input:
        vcf = config["results"] + "structural_variants/sniffles/LRS.bcf",
    output:
        bcf = config["results"] + "structural_variants/sniffles/LRS.{SV_type}.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools view {input.vcf} \
            -i "INFO/SVTYPE='{wildcards.SV_type}'" \
            -f \
            -Ob \
            -o {output.bcf} \
        """

## Mappability map

# Should try to merging this rule into the next
rule chop_ref_fna:
    """Chop reference FASTA into R1 and R2 FASTQ.

    Note that this rule will take up a lot of storage, about 206 times the size of the reference genome itself."""
    input:
        ref_fasta = config["ref_fasta"],
    output:
        R1 = temp(config["results"] + "structural_variants/delly/mappability/chopped_ref_fna.R1.fq.gz"),
        R2 = temp(config["results"] + "structural_variants/delly/mappability/chopped_ref_fna.R2.fq.gz"),
    params:
        R1 = config["results"] + "structural_variants/delly/mappability/chopped_ref_fna.R1",
        R2 = config["results"] + "structural_variants/delly/mappability/chopped_ref_fna.R2",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly2.yaml"
    shell: """
        dicey chop {input.ref_fasta} \
            --fq1 {params.R1} \
            --fq2 {params.R2} \
        """

rule align_chopped_ref_fna:
    """Align chopped reference FASTA."""
    input:
        ref_fasta = config["ref_fasta"],
        ref_indices = multiext(config["ref_fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
        R1 = config["results"] + "structural_variants/delly/mappability/chopped_ref_fna.R1.fq.gz",
        R2 = config["results"] + "structural_variants/delly/mappability/chopped_ref_fna.R2.fq.gz",
    output:
        bam = config["results"] + "structural_variants/delly/mappability/self_aligned.bam",
    params:
        bwa_threads = 42,
        samtools_threads = 8,
    threads: 50
    resources: nodes = 50
    conda: "../envs/common.yaml"
    shell: """
        bwa-mem2 {input.ref_fasta} {input.R1} {input.R2} \
            -t {params.bwa_threads} \
        | samtools sort \
            -@ {params.samtools_threads} \
            -o {output.bam} \
        """

rule create_mappability_map:
    """Create mappability map for calling CNVs."""
    input:
        bam = config["results"] + "structural_variants/delly/mappability/self_aligned.bam",
        bai = config["results"] + "structural_variants/delly/mappability/self_aligned.bam.bai",
    output:
        map = config["results"] + "structural_variants/delly/mappability/self_mapped.fa.gz",
    params:
        map = config["results"] + "structural_variants/delly/mappability/self_mapped.fa",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly2.yaml"
    # Not sure if --chromosome is actually necessary since input.bam will only be one chromosome anyway.
    shell: """
        dicey mappability2 {input.bam} \
            -o {output.map}; \
        gunzip {output.map} && bgzip {params.map} \
        """

## Copy number variation calling

rule call_CNVs:
    """Call CNVs per sample."""
    input:
        bam = config["results"] + "alignments/recalibrated/{sample}.bam",
        ref_fasta = config["ref_fasta"],
        map = config["results"] + "structural_variants/delly/mappability/map.fa.gz",
        fai = config["results"] + "structural_variants/delly/mappability/map.fa.gz.fai",
    output:
        bcf = config["results"] + "structural_variants/delly/CNVs/{sample}.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly2.yaml"
    shell: """
        delly cnv {input.bam} \
            -g {input.ref_fasta} \
            -l {input.map} \
            -o {output.bcf} \
        """

rule merge_CNVs:
    """Merge CNVs."""
    input:
        bcfs = expand(config["results"] + "structural_variants/delly/CNVs/{sample}.bcf", sample=SAMPLES),
    output:
        bcf = config["results"] + "structural_variants/delly/CNVs/{dataset}.ungenotyped.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly2.yaml"
    shell: """
        delly merge {input.bcfs} \
            --cnvmode \
            --pass \
            --minsize 1000 \
            --maxsize 100000 \
            -o {output.bcf} \
        """

rule genotype_merged_CNVs:
    """Genotype merged BCFs with copy number variants."""
    input:
        # bam = config["results"] + "alignments/recalibrated/{sample}.bam",
        bams = collect_runs_from_sample,  # From variant_calling.smk
        bams_idx = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_sample(wildcards))),
        bcf = config["results"] + "structural_variants/delly/merged/{dataset}.ungenotyped.bcf",
        ref_fasta = config["ref_fasta"],
        map = config["results"] + "structural_variants/delly/mappability/map.fa.gz",
        fai = config["results"] + "structural_variants/delly/mappability/map.fa.gz.fai",
    output:
        bcf = config["results"] + "structural_variants/delly/CNVs/all_sites/{dataset}.{sample}.genotyped.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly2.yaml"
    shell: """
        delly cnv <(samtools merge {input.bams} -o -) \
            --segmentation \
            -v {input.bcf} \
            -g {input.ref_fasta} \
            -m {input.map} \
            -o {output.bcf} \
        """

rule merge_all_sites_CNVs:
    """Merge samples with sites from all other samples."""
    input:
        bcfs = lambda wildcards: expand(
            config["results"] + "structural_variants/delly/CNVs/all_sites/{dataset}.{sample}.genotyped.bcf",
                dataset=wildcards.dataset,
                sample=SAMPLES),
    output:
        bcf = config["results"] + "structural_variants/delly/CNVs/merged/{dataset}.genotyped.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/common.yaml"
    shell: """
        bcftools merge {input.bcfs} \
            -m id \
            -Ob \
            -o {output.bcf} \
        """

rule filter_CNVs:
    input:
        bcf = config["results"] + "structural_variants/delly/CNVs/merged/{dataset}.genotyped.bcf",
        csi = config["results"] + "structural_variants/delly/CNVs/merged/{dataset}.genotyped.bcf.csi",
    output:
        bcf = config["results"] + "structural_variants/delly/CNVs/merged/{dataset}.genotyped.filtered.bcf",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly2.yaml"
    shell: """
        delly classify {input.bcf} \
            --filter germline \
            -o {output.bcf} \
        """

## Read-deapth profile

rule read_depth_profile:
    """Generate a read-depth profile of a sample."""
    input:
        bams = collect_runs_from_sample,  # From variant_calling.smk
        bams_idx = lambda wildcards: list(map(lambda bam: bam + ".bai", collect_runs_from_sample(wildcards))),
        map = config["results"] + "structural_variants/delly/mappability/self_mapped.fa.gz",
        ref_fasta = config["ref_fasta"],
        # Indexed ref_fasta
    output:
        read_depth = config["results"] + "structural_variants/delly/read_depth/{sample}.cov.gz",
    threads: 1
    resources: nodes = 1
    conda: "../envs/delly2.yaml"
    shell: """
        delly cnv {input.bams} \
            -a \
            -g {input.ref_fasta} \
            -m {input.map} \
            -o {output.read_depth} \
        """

## Comparing LRS and SRS

# As of Python 3.12, this rule runs, but gives a SyntaxWarning for line line: "-p '[A,C,T,G,N,\[,\]]' \"
# This warning appears every time snakemake is run, so this section is commented out.
# rule VCF_BND_to_BED:
#     """Compare LRS WGS to SRS WGS."""
#     wildcard_constraints:
#         seq = "LRS|WGS",
#         SV_type = "BND",
#     input:
#         vcf = "/master/abagwell/variant-analysis/results/rhesus/structural_variants/delly/per_sample/{seq}{sample_id}.vcf.gz",
#         #vcf = "/data/RHESUS/FASTQ/WGS/LRS/04.Result_X202SC23050126-Z01-F002.Macaca_mulatta_20231103/04.SV_VarDetect/LRS{sample_id}_r/LRS{sample_id}_r.sv.vcf.gz",
#     output:
#         bed = "/master/abagwell/workspace/LRS_comparison/{seq}{sample_id}.{SV_type}.bed",
#         #bed = "/master/abagwell/workspace/LRS_comparison/LRS{sample_id}.{SV_type}.bed",
#     conda: "../envs/common.yaml"
#     # Filter by SV type and "PASS"
#     # Remove metadata
#     # Remove "IMPRECISE" variants
#     # Keep only CHROM, POS, and ALT columns
#     # Strip out extra characters
#     # Merge columns (starting and ending breakpoint)
#     # Sort by chromosome and position
#     # Convert to BED (adds one to position, might technically be offset from actual values)
#     shell: """
#         bcftools view {input.vcf} \
#             -i "INFO/SVTYPE='{wildcards.SV_type}'" \
#             -f PASS \
#         | grep \
#             -v ^## \
#         | csvtk cut \
#             -t \
#             -C$ \
#             -f'#CHROM,POS,ALT' \
#         | csvtk replace \
#             -t \
#             -C$ \
#             -fALT \
#             -p '[A,C,T,G,N,\[,\]]' \
#             -r '' \
#         | sed '1d;s/\\t/\\n/2;s/:/\\t/1;' \
#         | csvtk sort \
#             -t \
#             -H \
#             -k 1:N,2:N \
#         | awk 'BEGIN {{OFS="\\t"}} {{print $1,$2,$2+11}}' \
#         > {output.bed} \
#         """

rule VCF_nonBND_to_BED:
    """Compare LRS WGS to SRS WGS."""
    wildcard_constraints:
        seq = "LRS|WGS",
        SV_type = "DEL|DUP|INS|INV",
    input:
        vcf = "/master/abagwell/variant-analysis/results/rhesus/structural_variants/delly/per_sample/{seq}{sample_id}.vcf.gz",
        #vcf = "/data/RHESUS/FASTQ/WGS/LRS/04.Result_X202SC23050126-Z01-F002.Macaca_mulatta_20231103/04.SV_VarDetect/LRS{sample_id}_r/LRS{sample_id}_r.sv.vcf.gz",
    output:
        bed = "/master/abagwell/workspace/LRS_comparison/{seq}{sample_id}.{SV_type}.bed",
        #bed = "/master/abagwell/workspace/LRS_comparison/LRS{sample_id}.{SV_type}.bed",
    conda: "../envs/common.yaml"
    # Filter by SV type and "PASS"
    # Remove metadata
    # Remove "IMPRECISE" variants
    # Keep only CHROM, POS, and ALT columns
    # Strip out extra characters
    # Merge columns (starting and ending breakpoint)
    # Sort by chromosome and position
    # Convert to BED (adds one to position, might technically be offset from actual values)
    shell: """
        bcftools view {input.vcf} \
            -i "INFO/SVTYPE='{wildcards.SV_type}'" \
            -f PASS \
        | grep \
            -v ^## \
        | csvtk cut \
            -t \
            -C$ \
            -f'#CHROM,POS' \
        | sed '1d;' \
        | csvtk sort \
            -t \
            -H \
            -k 1:N,2:N \
        | awk 'BEGIN {{OFS="\\t"}} {{print $1,$2,$2+11}}' \
        > {output.bed} \
        """
        
rule intersect_LRS_to_SRS:
    """Compare LRS WGS to SRS WGS."""
    input:
        LRS = "/master/abagwell/workspace/LRS_comparison/LRS{sample_id}.{SV_type}.bed",
        WGS = "/master/abagwell/workspace/LRS_comparison/WGS{sample_id}.{SV_type}.bed",
    output:
        LRS = "/master/abagwell/workspace/LRS_comparison/{sample_id}.{SV_type}.loj_LRS.tsv",
        WGS = "/master/abagwell/workspace/LRS_comparison/{sample_id}.{SV_type}.loj_WGS.tsv",
    shell: """
        bedtools intersect \
            -a {input.LRS} \
            -b {input.WGS} \
            -loj \
        > {output.LRS}; \
        bedtools intersect \
            -a {input.WGS} \
            -b {input.LRS} \
            -loj \
        > {output.WGS}; \
        """

rule compile_intersection_stats:
    """Compile intersection stats."""
    input:
        LRS = "/master/abagwell/workspace/LRS_comparison/{sample_id}.{SV_type}.loj_LRS.tsv",
        WGS = "/master/abagwell/workspace/LRS_comparison/{sample_id}.{SV_type}.loj_WGS.tsv",
    output:
        stats = "/master/abagwell/workspace/LRS_comparison/{sample_id}.{SV_type}.stats",
    shell: """
        echo "ID\tSV_Type\tLRS_Total\tMatches\tWGS_Total" >> {output.stats}; \

        LRS_TOTAL=$(wc -l {input.LRS} | cut -d ' ' -f 1); \
        MATCHES=$(cat {input.LRS} | grep -e "-1" -v | wc -l | cut -d ' ' -f 1); \
        WGS_TOTAL=$(wc -l {input.WGS} | cut -d ' ' -f 1); \

        echo "{wildcards.sample_id}\t{wildcards.SV_type}\t$LRS_TOTAL\t$MATCHES\t$WGS_TOTAL" >> {output.stats}; \
        """

rule concat_intersection_stats:
    input:
        stats = expand("/master/abagwell/workspace/LRS_comparison/{sample_id}.{SV_type}.stats", sample_id=["31236", "32351", "32592", "38454"], SV_type=["BND", "DEL", "DUP", "INS", "INV"]),
    output:
        merged = "/master/abagwell/workspace/LRS_comparison/merged.stats",
    shell: """
        csvtk concat {input.stats} \
            -t \
            -C$ \
        | csvtk sort \
            -t \
            -k SV_Type,ID \
        > {output.merged} \
        """

