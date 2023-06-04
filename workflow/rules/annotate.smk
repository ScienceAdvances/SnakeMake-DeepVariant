rule annotate_variants:
    input:
        calls="results/called/all.vcf.gz",  # .vcf, .vcf.gz or .bcf
        cache=genome_prefix+"vep_cache",  # can be omitted if fasta and gff are specified
        plugins=genome_prefix+"vep_plugins",
        # optionally add reference genome fasta
        # fasta=gatk_dir+"genome.fa",
        # fai=gatk_dir+"genome.fa.fai", # fasta index
        # gff=gatk_dir+"annotation.gff.gz",
        # csi=gatk_dir+"annotation.gff.gz.csi", # tabix index
        # add mandatory aux-files required by some plugins if not present in the VEP plugin directory specified above.
        # aux files must be defined as following: "<plugin> = /path/to/file" where plugin must be in lowercase
        # revel = path/to/revel_scores.tsv.gz
    output:
        calls="results/vep/all.vcf.gz",
        stats=report("report/vep_report.html",caption="../report/Annotation.rst",category="Annotation")
    params:
        plugins=config["vep"]["plugins"],
        extra=config["vep"]["extra"],
    log:
        "logs/vep/annotate.log"
    threads: 16
    wrapper:
        config["warpper_mirror"]+"bio/vep/annotate"