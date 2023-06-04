rule get_genome:
    output:
        genome_prefix+"genome.fa"
    log:
        "logs/ref/get_genome.log"
    retries: 50
    params:
        species=config["genome"].get("species"),
        datatype=config["genome"].get("datatype"),
        build=config["genome"].get("build"),
        release=config["genome"].get("release"),
    cache: True
    wrapper:
        config["warpper_mirror"]+"bio/reference/ensembl-sequence"

rule genome_faidx:
    input:
        genome_prefix+"genome.fa"
    output:
        genome_prefix+"genome.fa.fai"
    cache: True
    log:
        "logs/ref/genome_faidx.log"
    wrapper:
        config["warpper_mirror"]+"bio/samtools/faidx"

rule picard_create_dict:
    input:
        genome_prefix+"genome.fa"
    output:
        genome_prefix+"genome.dict"
    log:
        "logs/ref/picard_create_dict.log"
    resources:
        mem_mb=1024
    cache: True
    wrapper:
        config["warpper_mirror"]+"bio/picard/createsequencedictionary"

rule bwa_mem2_index:
    input:
        genome_prefix+"genome.fa"
    output:
        multiext(
            genome_prefix+"genome.fa",
            ".0123",
            ".amb",
            ".ann",
            ".bwt.2bit.64",
            ".pac",
        ),
    log:
        "logs/ref/bwa_mem2_index.log"
    params:
        bwa="bwa-mem2"
    threads: 16
    cache: True
    wrapper:
        config["warpper_mirror"]+"bio/bwa-memx/index"

rule get_vep_cache:
    output:
        directory(genome_prefix+"vep_cache"),
    params:
        species=config["genome"].get("species"),
        build=config["genome"].get("build"),
        release=config["genome"].get("release"),
    retries: 50
    log:
        "logs/ref/get_vep_cache.log"
    cache: True
    wrapper:
        config["warpper_mirror"]+"bio/vep/cache"

rule get_vep_plugins:
    output:
        directory(genome_prefix+"vep_plugins")
    retries: 50
    params:
        config["genome"].get("release")
    log:
        "logs/ref/get_vep_plugins.log"
    cache: True
    shell:
        "git clone -b release/{params[0]} --single-branch https://jihulab.com/BioQuest/vep_plugins.git {output[0]} 2> {log}"