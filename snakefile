with open('SRR_Acc_List.txt') as f:
    accession = [line.rstrip() for line in f]

accession = [i for i in accession if not '#' in i]

configfile: "config.yml"
containerized: "container.sif"

rule all:
    input:
        expand("data/fastq/{sra}_{direction}.fastq.gz", sra = accession, direction = ["1","2"]),
        expand("results/fastqc/{sample}_{direction}_fastqc.html", sample = accession, direction = ["1", "2"]),
        expand("results/fastqc/{sample}_{direction}_processed_fastqc.html", sample = accession, direction = ["1", "2"]),
        expand("processed_data/fastq/{sample}_{direction}_processed.fastq", sample = accession, direction = ["1", "2"]),
        expand("results/aligned_reads/{sample}_{genome}.bam.bai", sample = accession, genome = ['allgenomes_sorted', 'cmv_extract', 'human_extract']),
        expand("results/aligned_reads/{sample}_{genome}_extract.bam", sample = accession, genome = ['cmv', 'spikein', 'human']),
        expand("results/bed/{sample}_{genome}.bed", sample = accession, genome = ['cmv', "human"]),
        expand("results/tracks/{sample}_{genome}_{direction}.bg", sample = accession, genome = ['cmv', "human"], direction = ['for', 'rev']),
        "results/multiqc/multiqc_report.html"

include: "rules/get_fastq.smk"
include: "rules/process_fastq.smk"
include: "rules/quality_control.smk"
include: "rules/align_reads.smk"
include: "rules/create_bed.smk"
include: "rules/create_track.smk"


