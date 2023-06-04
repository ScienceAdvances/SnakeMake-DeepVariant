DeepVariant Pipeline summary

SbankeMake Workflow for DeepVariant calling

=============================================
Reference
=============================================

Genome refence for {{ snakemake.config["genome"]["species"] }} build {{ snakemake.config["genome"]["build"] }} release {{ snakemake.config["genome"]["release"] }} was download from Ensembl_

=============================================
Align
=============================================

1. Adapter trimming (fastp_)
2. Aligner (`BWA mem2`_)
3. Mark duplicates (samblaster_)
4. Merge CRAMs of every sample, repesectly (Picard_)
5. Create CRAM index (samtools_)

=============================================
Quality control report
=============================================

1. Fastp report (MultiQC_)
2. Alignment report (MultiQC_)

=============================================
Germline variants
=============================================

1. Predict variants (DeepVariant_)
2. Joint variant calling (GLnexus_)

=============================================
Annotation
=============================================

Annotate variant calls with VEP (VEP_)

.. _Ensembl: https://asia.ensembl.org/index.html
.. _VEP: https://www.ensembl.org/info/docs/tools/vep/index.html
.. _fastp: https://github.com/OpenGene/fastp
.. _BWA mem2: http://bio-bwa.sourceforge.net
.. _samblaster: https://github.com/GregoryFaust/samblaster
.. _MultiQC: https://multiqc.info
.. _samtools: http://www.htslib.org
.. _GLnexus: https://github.com/dnanexus-rnd/GLnexus
.. _Picard: https://broadinstitute.github.io/picard
.. _DeepVariant: https://github.com/google/deepvariant