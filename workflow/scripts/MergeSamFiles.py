"""Snakemake wrapper for picard MergeSamFiles."""

__author__ = "Victor Wang"
__copyright__ = "Copyright 2023, Victor Wang"
__email__ = "victor@bioquest.cn"
__license__ = "Apache License 2.0"


import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

picard_extra = snakemake.params.get("picard_extra", "")
samtools_extra = snakemake.params.get("samtools_extra", "")
java_opts = get_java_opts(snakemake)
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
maps=snakemake.input.get("maps",False)
if not maps:
    raise ValueError("input file is None")
if isinstance(maps,list) & len(maps)>1:
    maps = " ".join("--INPUT {}".format(x) for x in maps)
    run_merger = True
else:
    run_merger = False

ref=snakemake.input.get("ref","")
if ref:
    ref=f"--REFERENCE_SEQUENCE {ref}"
create_idx=snakemake.output.get("idx","")
if create_idx:
    create_idx="--CREATE_INDEX true"
threads=snakemake.threads
threads = "" if threads <= 1 else " -@ {} ".format(snakemake.threads - 1)
if threads:
    use_threads="--USE_THREADING true"

if run_merger:
    with tempfile.TemporaryDirectory() as tmpdir:
        shell(
            "picard MergeSamFiles"
            " {java_opts} {extra}"
            " {ref}"
            " {maps}"
            " --TMP_DIR {tmpdir}"
            " --OUTPUT {snakemake.output.map}"
            " {create_idx}"
            " --SORT_ORDER coordinate"
            " {use_threads}"
            " {log}"
        )
else:
    shell(
        "cp {snakemake.input.maps} {snakemake.output.map}" 
    )
    shell(
        "samtools index {threads} {samtools_extra} {snakemake.input.maps} {snakemake.output.idx}"
    )
