rule DeepVariant:
    input:
        reads="results/aligned/{s}.merged.cram",
        ref=genome_prefix+"genome.fa",
        fai=genome_prefix+"genome.fa.fai",
        intervals=config["fastqs"]["restrict_regions"],
    output:
        vcf="results/called/{s}.vcf.gz",
        gvcf="results/called/{s}.g.vcf.gz",
        intermediate_results_dir=directory("results/called/{s}.intermediate_results")
    params:
        omics=config["fastqs"]["omics"],
    threads: 32
    container: 
        "docker://hub-mirror.c.163.com/google/deepvariant:1.5.0"
    log:
        "logs/DeepVariant/{s}.log"
    script:
        "../scripts/DeepVariant.py"

rule GLnexus:
    input:
        gvcfs=expand("results/called/{s}.g.vcf.gz",s=samples.Sample),
        intervals=config["fastqs"]["restrict_regions"],
    output:
        "results/called/all.vcf.gz"
    params:
        omics=config["fastqs"]["omics"],
    threads: 32
    log:
        "logs/GLnexus.log",
    conda:
        "../envs/GLnexus.yaml"
    # container:
    #     "docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.3"
    script:
        "../scripts/GLnexus.py"