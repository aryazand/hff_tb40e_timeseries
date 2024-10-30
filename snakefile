with open('SRR_Acc_List.txt') as f:
    accession = [line.rstrip() for line in f]
snakemake
accession = [i for i in accession if not '#' in i]

configfile: "config.yml"
containerized: "container.sif"

rule all:
    input:
        expand("processed_data/fastq/{sample}_{direction}_processed.fastq.gz", sample = accession, direction = ["1", "2"]),
        expand("results/aligned_reads/{sample}_{genome}.bam.bai", sample = accession, genome = ['allgenomes_sorted', 'cmv_extract', 'human_extract']),
        expand("results/aligned_reads/{sample}_{genome}_extract.bam", sample = accession, genome = ['cmv', 'spikein', 'human']),
        expand("results/bed/{sample}_{genome}.bed", sample = accession, genome = ['cmv', "human"]),
        expand("results/tracks/{sample}_{genome}_{direction}.{track_type}", sample = accession, genome = ['cmv', "human"], direction = ['for', 'rev'], track_type = config['create_track']['track_type']),
        "results/QC/multiqc/multiqc_report.html"

include: "rules/get_fastq.smk"
include: "rules/process_fastq.smk"
include: "rules/quality_control.smk"
include: "rules/align_reads.smk"
include: "rules/create_bed.smk"
include: "rules/create_track.smk"
