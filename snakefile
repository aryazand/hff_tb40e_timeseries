#########################
# Define Environment
#########################

import os
import peppy
import pandas as pd

configfile: "config.yml"
# containerized: "container.sif"

include: "rules/common.smk"

#############################
# Define directories
#############################

DATA_DIR = "data"
RAWFASTQ_DIR = os.path.join(DATA_DIR, "fastq")
GENOMES_DIR = os.path.join(DATA_DIR, "genomes")
RESULTS_DIR = config["results"]
TRIMMED_DIR = os.path.join(RESULTS_DIR, "trimmed")
DEDUPPED_DIR = os.path.join(RESULTS_DIR, "dedupped")
ALIGNMENT_DIR = os.path.join(RESULTS_DIR, "alignments")
LOGS_DIR = os.path.join(RESULTS_DIR, "logs")

QC_DIR = os.path.join(RESULTS_DIR, "qc")
FASTQC_RAW_DIR = os.path.join(QC_DIR, "fastqc_raw")
FASTQC_PROCESSED_DIR = os.path.join(QC_DIR, "fastqc_processed")
QC_DIR_TRIMMING = os.path.join(QC_DIR, "trimming_reports")
QC_DIR_ALIGNMENT = os.path.join(QC_DIR, "alignment_reports")
QC_DIR_DEDUP = os.path.join(QC_DIR, "dedup")
QC_DIR_UMIEXTRACT = os.path.join(QC_DIR, "umi_extraction")
MULTIQC_DIR = os.path.join(QC_DIR, "multiqc")

#############################
# Load samples and metadata
#############################

# Load project PEP
project = peppy.Project("data/sample_metadata/PEP.yaml")

# Get sample metadata
sample_table = project.sample_table

# Only keep samples with flavo
sample_table = sample_table[sample_table['sample_name'].str.contains('flavo')]

# For testing, only keep one sample
#sample_table = sample_table[sample_table['srr'] == "SRR16301852"]
#sample_table["srr"] = "SRR5660034"	

# Create a genomes key
dict_list = []
for i in range(len(sample_table)):
    genome_names = sample_table['genome_names'].iloc[i].split(", ")
    genome_accessions = sample_table['genome_accessions'].iloc[i].split(", ")
    dict_list.append(dict(zip(genome_names, genome_accessions)))

for i in range(len(dict_list)-1):
    dict_list[i].update(dict_list[i+1])

GENOMES = dict_list[0]

#############################
# Create target files 
#############################

# Uses snakemake expand and zips the sample names with the genome names
# to create a list of bam files
all_genomes_bam = expand(os.path.join(ALIGNMENT_DIR,"{sample}_allgenomes.bam"), sample = sample_table['sample_name'])

######################
# Define output files
#####################

rule all:
    input:
        expand(os.path.join(ALIGNMENT_DIR,"{sample}_allgenomes.dedup.bam"), sample = sample_table['sample_name']),
        expand(os.path.join(ALIGNMENT_DIR,"{sample}_allgenomes.dedup.bam.bai"), sample = sample_table['sample_name']),
        expand(os.path.join(ALIGNMENT_DIR,"{sample}_{genome}.bam"), sample = sample_table['sample_name'], genome = GENOMES.keys()),
        expand(os.path.join(ALIGNMENT_DIR,"{sample}_{genome}.bam.bai"), sample = sample_table['sample_name'], genome = GENOMES.keys())

rule multiqc:
    input:
        os.path.join(MULTIQC_DIR, "multiqc_report.html")

include: "rules/get_fastq.smk"
include: "rules/process_fastq.smk"
include: "rules/quality_control.smk"
include: "rules/align_reads.smk"
#include: "rules/create_bed.smk"
include: "rules/get_genomic_data.smk"
