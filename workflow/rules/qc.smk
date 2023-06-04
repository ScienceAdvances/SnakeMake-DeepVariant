rule fastp_multiqc:
    input:
        expand("report/{s}{u}.fastp.json", s=samples.Sample,u=samples.Unit)
    output:
        report(
            "report/fastp_multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality control"
        ),
    log:
        "logs/qc/fastp_multiqc.log"
    wrapper:
        config["warpper_mirror"]+"bio/multiqc"

rule samtools_stats:
    input:
        "results/aligned/{s}.cram",
        bed=config["fastqs"]["restrict_regions"],#Optional input, specify target regions
    output:
        temp("results/aligned/{s}.txt")
    log:
        "logs/qc/{s}.samtools_stats.log"
    wrapper:
        config["warpper_mirror"]+"bio/samtools/stats"

rule align_multiqc:
    input:
        expand(["results/aligned/{s}.txt"],s=samples.Sample),
    output:
        report(
            "report/align_multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality control",
        ),
    log:
        "logs/qc/align_multiqc.log",
    wrapper:
        config["warpper_mirror"]+"bio/multiqc"