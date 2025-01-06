#########################
# Define Environment
#########################

import os
import pandas as pd
import pathlib
from snakemake.utils import Paramspace

configfile: "config.yml"
containerized: "container.sif"

include: "rules/common.smk"

#########################
# Define global variables
#########################

samples = Paramspace(pd.read_csv("sample_metadata.csv"))
sample_names = set(samples.sample_name)
sample_ext = get_fastq_extention(samples)

SPECIES = list(config["genomes"].keys())
SPIKEIN_SPECIES = list(config['spikein'].keys())
GENOMES = {}
ACCESSION = {}

for specie in SPECIES:
    GENOMES[specie] = config["genomes"][specie]["genome_name"]
    ACCESSION[specie] = config["genomes"][specie]["accession"]

for specie in SPIKEIN_SPECIES:
    GENOMES[specie] = config["spikein"][specie]["genome_name"] 
    ACCESSION[specie] = config["spikein"][specie]["accession"] 

BIGWIGS_FOR_USCSC = []
USCSC_HUB = []

for species in SPECIES:
    for sample in sample_names:
        for direction in ['for', 'rev']:
            genome = config["genomes"][species]["genome_name"]
            bw = f"results/UCSCGenomeBrowser/{species}/{genome}/bw/{sample}_{species}_{direction}.bw"
            BIGWIGS_FOR_USCSC.append(bw)

for species in SPECIES:
        genome = config["genomes"][species]["genome_name"]
        hub = f"results/UCSCGenomeBrowser/{species}/{genome}/trackDb.txt"
        USCSC_HUB.append(hub)

######################
# Define output files
#####################

rule all:
    input:
        expand("results/aligned_reads/{sample}_{genome}.bam.bai", 
            sample = sample_names, 
            genome = ['allgenomes_sorted', 'cmv_extract', 'human_extract']),
        expand("results/aligned_reads/{sample}_{genome}_extract.bam", 
            sample = sample_names, 
            genome = SPECIES + SPIKEIN_SPECIES),
        expand("results/bed/{sample}_{genome}.bed", 
            sample = sample_names, 
            genome = SPECIES),
        expand("results/tracks/{sample}_{genome}_{direction}.{track_type}", 
            sample = sample_names, 
            genome = SPECIES, 
            direction = ['for', 'rev'], 
            track_type = config['create_track']['track_type']),
        expand("results/tracks/{sample}_{genome}_{direction}.fiveprime.bw", 
            sample = sample_names, 
            genome = SPECIES, 
            direction = ['for', 'rev']),
        "results/QC/multiqc/multiqc_report.html",
        BIGWIGS_FOR_USCSC,
        expand("results/UCSCGenomeBrowser/{species}/genomes.txt", 
            species = SPECIES),
        expand("results/UCSCGenomeBrowser/{species}/hub.txt", 
            species = SPECIES),
        USCSC_HUB

include: "rules/get_fastq.smk"
include: "rules/process_fastq.smk"
include: "rules/quality_control.smk"
include: "rules/align_reads.smk"
include: "rules/create_bed.smk"
include: "rules/create_track.smk"
include: "rules/create_ucsc_hub.smk"
include: "rules/get_genomic_data.smk"