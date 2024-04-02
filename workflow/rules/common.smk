import numpy as np
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("7.25.0")

report: "../report/workflow.rst"

#=====================================================
# validate config.yaml file and samples.csv file
#=====================================================

configfile: "config/config.yaml"

validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], dtype=str,sep='\t',header=0).fillna(value="")

if not "Unit" in samples.columns:
    samples.loc[:,"Unit"]=""
samples.loc[:,"Unit"] = [f".{x}" if x else x for x in samples.Unit]
samples.set_index(keys=["Sample", "Unit"], drop=False,inplace=True)

samples.index = samples.index.set_levels(
    [i.astype(str) for i in samples.index.levels]
)  # enforce str in index

validate(samples, schema="../schemas/samples.schema.yaml")

# if units are not none, add a . prefix
fastqs = config["fastqs"].get("dir")
if config["fastqs"].get("pe"):
    fq1=[f"{fastqs}/{x}{y}.1.fastq.gz" for x,y in zip(samples.Sample,samples.Unit)]
    fq2=[f"{fastqs}/{x}{y}.2.fastq.gz" for x,y in zip(samples.Sample,samples.Unit)]
    samples.insert(loc=0,column="fq2",value=fq2)
    samples.insert(loc=0,column="fq1",value=fq1)
else:
    fq1=[f"{fastqs}/{x}{y}.fastq.gz" for x,y in zip(samples.Sample,samples.Unit)]
    samples.insert(loc=0,column="fq1",value=fq1)

#=====================================================
# Helper functions
#=====================================================

def get_genome_prefix(config=config):
    """genome files prefix"""
    g=config["genome"]
    p="{0}/{1}_{2}_{3}_".format(g['dir'],g['species'].capitalize(),g['build'],g['release'])
    return p
genome_prefix=get_genome_prefix(config)

def get_dedup_extra(config=config):
    if config["fastqs"]["pe"]:
        return "-M"
    else:
        return "-M --ignoreUnmated"

# def get_captured_regions(config=config):
#     c=pd.read_csv(config["fastqs"]["restrict_regions"],sep='\t',header=None,dtype=str)
#     return [f"{x}:{y}-{z}" if y else x for _,x,y,z in c.itertuples()]

def get_fastq(wildcards):
    """Get fastq files of given sample and unit."""
    fastqs = samples.loc[(wildcards.s, wildcards.u), ]
    if config["fastqs"].get("pair_end"):
        return [fastqs.fq1, fastqs.fq2]
    return [fastqs.fq1]

def get_trimmed_fastq(wildcards):
    """Get trimmed reads of given sample and unit."""
    if config["fastqs"].get("pe"):
        # paired-end sample
        return expand("results/trimmed/{{s}}{{u}}.{_}.fastq.gz",_=["1","2"])
    # single end sample
    return ["results/trimmed/{s}{u}.fastq.gz"]

def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{s}{u}\tSM:{s}\tLB:{s}\tPL:{pl}' -M".format(
        s=wildcards.s,
        u=wildcards.u,
        pl=config["fastqs"].get("platform")
    )

def get_cram(wildcards):
    """Get all aligned reads of given sample."""
    u=samples.loc[wildcards.s,"Unit"]
    if u.all(axis=None):
        return expand(
            ["results/aligned/{{s}}{u}.cram"],
            u=u
        )
    else:
        return "results/aligned/{s}.cram"
#=====================================================
# Wildcard constraints
#=====================================================

wildcard_constraints:
    s="|".join(samples.index.get_level_values(0)),
    u="|".join(samples.index.get_level_values(1)),