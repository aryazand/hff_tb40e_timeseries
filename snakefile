#########################
# Define Environment
#########################

import os
import peppy
from snakemake.utils import Paramspace

configfile: "config.yml"
# containerized: "container.sif"

include: "rules/common.smk"

#############################
# Define directories
#############################

RESULTS_DIR = config["results"]
TRIMMED_DIR = os.path.join(RESULTS_DIR, "trimmed")
DEDUPPED_DIR = os.path.join(RESULTS_DIR, "dedupped")
ALIGNMENT_DIR = os.path.join(RESULTS_DIR, "alignments")
LOGS_DIR = os.path.join(RESULTS_DIR, "logs")
QC_DIR = os.path.join(RESULTS_DIR, "qc")

#############################
# Load samples and metadata
#############################

# Load project PEP
project = peppy.Project("data/sample_metadata/PEP.yaml")

# Get sample metadata
sample_table = project.sample_table

# Only keep samples with flavo
sample_table = sample_table[sample_table['sample_name'].str.contains('flavo')]
sample_table = sample_table[sample_table['srr'] == "SRR16301852"]
sample_table = Paramspace(sample_table)

#############################
# Create target files 
#############################

# Uses snakemake expand and zips the sample names with the genome names
# to create a list of bam files
all_genomes_bam = expand("results/aligned_reads/{sample}_allgenomes.bam", sample = sample_table['sample_name'])

######################
# Define output files
#####################

rule all:
    input:
        expand(os.path.join(TRIMMED_DIR, "{sample}_{direction}_val_{direction}.fq.gz"), 
            sample = sample_table.sample_name, 
            direction = ['1', '2']),
#         expand("results/aligned_reads/{sample}_{genome}.bam.bai", 
#             sample = sample_names, 
#             genome = ['allgenomes_sorted', 'cmv_extract', 'human_extract']),
#         expand("results/aligned_reads/{sample}_{genome}_extract.bam", 
#             sample = sample_names, 
#             genome = SPECIES + SPIKEIN_SPECIES),
#         expand("results/bed/{sample}_{genome}.bed", 
#             sample = sample_names, 
#             genome = SPECIES),
#         "results/QC/multiqc/multiqc_report.html"

include: "rules/get_fastq.smk"
include: "rules/process_fastq.smk"
# include: "rules/quality_control.smk"
# include: "rules/align_reads.smk"
# include: "rules/create_bed.smk"
# include: "rules/get_genomic_data.smk"
