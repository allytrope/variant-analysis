
env_path = "../../envs"

rule paragraph_idxdepth:
    """Find the average depth across genome for the manifest file for paraGRAPH."""
    input:
        alignment = config["results"] + "alignments/pbmm2/{batch}/{seq}{indiv}_{library}.bam",
        index = config["results"] + "alignments/pbmm2/{batch}/{seq}{indiv}_{library}.bam.bai",
        ref_fasta = "/master/abagwell/variant-analysis/resources/rhesus/ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa",
    output:
        idxdepth = config["results"] + "structural_variants/paraGRAPH/{batch}/{seq}{indiv}_{library}.txt"
    threads: 4  # Default is 48
    resources:
        nodes = 1
    conda: f"{env_path}/paragraph.yaml"
    shell: """
        idxdepth \
            --bam {input.alignment} \
            --reference {ref_fasta} \
            --threads {threads} \
        > {output.idxdepth} \
        """

# TODO: Edit based on how idxdepth output looks
rule paragraph_manifest:
    """Create the "manifest" file for paraGRAPH."""
    input:
        idxdepth = config["results"] + "structural_variants/paraGRAPH/{batch}/{seq}{indiv}_{library}.txt"
    output:
        tsv = config["results"] + "structural_variants/paraGRAPH/{dataset}.manifest.tsv"
    threads: 1
    resources:
        nodes = 1
    conda: "{env_path}/paragraph.yaml"
    shell: """
        
        """


rule paragraph:
    """Call structual variants on short-read data but using
    structural variants from long-read data as a reference."""
    input:
        candidate_vcf = config["results"] + "structural_variants/SURVIVOR/{dataset}.chr{chr}.vcf.gz",
        samples_file = config["results"] + "structural_variants/paraGRAPH/{dataset}.chr{chr}.samples.txt",
        ref_fasta = "/master/abagwell/variant-analysis/resources/rhesus/ref_fna/Macaca_mulatta.Mmul_10.dna.toplevel.fa",
    output:
        config["results"] + "structural_variants/paraGRAPH/genotypes.vcf.gz",
        config["results"] + "structural_variants/paraGRAPH/genotypes.json.gz",
        config["results"] + "structural_variants/paraGRAPH/variants.vcf.gz",
        config["results"] + "structural_variants/paraGRAPH/variants.json.gz",
    log:
        config["results"] + "structural_variants/paraGRAPH/grmpy.log"
    params:
        output_dir = config["results"] + "structural_variants/paraGRAPH",
    # Apparently paragraph actually runs on the python script "multigrmpy.py"
    # Actually, elsewhere though, there is a reference to "paragraph" and "grmpy" commands
    threads: 1
    resources:
        nodes = 1
    conda: f"{env_path}/paragraph.yaml"
    shell: """
        python /master/abagwell/tools/paraGRAPH/multigrmpy.py \
            -i {input.candidate_vcf} \
            -m {input.samples_file} \
            -r {input.ref_fasta} \
            -o {params.output_dir} \
        """

# rule delly: