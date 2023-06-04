__author__ = "Tetsuro Hisayoshi"
__copyright__ = "Copyright 2020, Tetsuro Hisayoshi"
__email__ = "hisayoshi0530@gmail.com"
__license__ = "MIT"

import os
import tempfile
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

# log_dir = os.path.dirname(snakemake.log[0])

intermediate_dir=snakemake.output.get("intermediate_results_dir", "")
if intermediate_dir:
    intermediate_results_dir=f"--intermediate_results_dir={intermediate_dir}"
else:
    intermediate_results_dir="--intermediate_results_dir={tmp_dir}"

intervals=snakemake.input.get("intervals", "")
if intervals:
    intervals=f"--regions={intervals}"

output_gvcf = snakemake.output.get("gvcf", "")
if output_gvcf:
    output_gvcf = f"--output_gvcf {output_gvcf}"

with tempfile.TemporaryDirectory() as tmp_dir:
    shell(
        "run_deepvariant"
        " --model_type={snakemake.params.omics}"
        " --ref={snakemake.input.ref}"
        # " --fai={snakemake.input.fai}"
        " --reads={snakemake.input.reads}"
        " {intervals}"
        " --output_vcf={snakemake.output.vcf}"
        " {output_gvcf}"
        " {intermediate_results_dir}"
        " --num_shards={snakemake.threads}"
        # " --logging_dir {log_dir}"
        " {extra} {log}"
    )