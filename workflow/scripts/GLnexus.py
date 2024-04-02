__author__ = "Victor Wang"
__copyright__ = ("Copyright 2023, Victor Wang")
__email__ = "victor@bioquest.cn"
__license__ = "Apache License 2.0"


import tempfile
from snakemake.shell import shell

# Extract arguments.
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

intervals = snakemake.input.get("intervals","")
if intervals:
    intervals=f"--bed {intervals}"

with tempfile.TemporaryDirectory() as tmp_dir:
    shell(
        "glnexus_cli"
        " --config {snakemake.params.config}"
        " {intervals}"
        " --dir {tmp_dir}glnexus_cli"
        " --threads {snakemake.threads}"
        " {snakemake.input.gvcfs}"
        " | bcftools view -"
        " | bgzip -@ {snakemake.threads} -c"
        " > {snakemake.output[0]} {log}"
    )