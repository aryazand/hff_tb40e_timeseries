include: "rules/common.smk"

import os
import pandas as pd
import pathlib
from snakemake.utils import Paramspace
samples = Paramspace(pd.read_csv("sample_metadata.csv"))
sample_names = set(samples.sample_name)
sample_ext = get_fastq_extention(samples)

configfile: "config.yml"
containerized: "container.sif"

rule all:
    input:
        expand("processed_data/fastq/{sample}_R{direction}_processed.fastq.gz", sample = sample_names, direction = ["1", "2"]),
        expand("results/aligned_reads/{sample}_{genome}.bam.bai", sample = sample_names, genome = ['allgenomes_sorted', 'cmv_extract', 'human_extract']),
        expand("results/aligned_reads/{sample}_{genome}_extract.bam", sample = sample_names, genome = ['cmv', 'spikein', 'human']),
        expand("results/bed/{sample}_{genome}.bed", sample = sample_names, genome = ['cmv', "human"]),
        expand("results/tracks/{sample}_{genome}_{direction}.{track_type}", sample = sample_names, genome = ['cmv', "human"], direction = ['for', 'rev'], track_type = config['create_track']['track_type']),
        expand("results/tracks/{sample}_{genome}_{direction}_fiverpime.bw", sample = sample_names, genome = ['cmv'], direction = ['for', 'rev']),
        "results/QC/multiqc/multiqc_report.html"

include: "rules/get_fastq.smk"
include: "rules/process_fastq.smk"
include: "rules/quality_control.smk"
include: "rules/align_reads.smk"
include: "rules/create_bed.smk"
include: "rules/create_track.smk"