rule align:
    input:
        reads=get_trimmed_fastq,
        reference=genome_prefix+"genome.fa",
        idx=multiext(genome_prefix+"genome.fa", ".0123", ".amb", ".bwt.2bit.64", ".ann",".pac"),
    output:
        "results/aligned/{s}{u}.cram" # Output can be .cram, .bam, or .sam
    log:
        "logs/align/{s}{u}.bwa-mem2.log"
    params:
        bwa="bwa-mem2", # Can be 'bwa-mem, bwa-mem2 or bwa-meme.
        extra=get_read_group,
        sort="picard",
        sort_order="coordinate",
        dedup=config['fastqs']['duplicates'], # Can be 'none' (default), 'mark' or 'remove'.
        dedup_extra=get_dedup_extra(),
        exceed_thread_limit=True,
        embed_ref=False,
    threads: 32
    wrapper:
        config["warpper_mirror"]+"bio/bwa-memx/mem"

rule MergeSamFiles:
    input:
        maps=get_cram,
        ref=genome_prefix+"genome.fa",
    output:
        map="results/aligned/{s}.merged.cram",
        idx="results/aligned/{s}.merged.cram.crai"
    threads: 32
    params:
        samtools_extra = "-c"
    log:
        "logs/align/MergeSamFiles/{s}.MergeSamFiles.log"
    conda: "../envs/MergeSamFiles.yaml"
    script:
        "../scripts/MergeSamFiles.py"