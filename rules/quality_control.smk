rule fastqc:
    input: 
        "data/fastq/{sample}.fastq"
    output: 
        "results/fastqc/{sample}_fastqc.html"
    conda:
        "../envs/proseq-qc.yml"
    threads: 5
    log:
        out = "log/fastqc_{sample}.out",
        err = "log/fastqc_{sample}.err"
    shell: 
        "fastqc {input} -t {threads} -o results/fastqc 2> {log.err} 1> {log.out}"


rule fastqc_on_processed_fastq:
    input: 
        "processed_data/fastq/{sample}_{direction}_processed.fastq.gz"
    output: 
        "results/fastqc/{sample}_{direction}_processed_fastqc.html"
    conda:
        "../envs/proseq-qc.yml"
    threads: 5
    log:
        out = "log/fastqc_{sample}_{direction}_processed.out",
        err = "log/fastqc_{sample}_{direction}_processed.err"
    shell: 
        "fastqc {input} -t {threads} -o results/fastqc 2> {log.err} 1> {log.out}"

rule run_mutliqc:
    input:
        expand("results/aligned_reads/{sample}_{genome}_extract.bam", sample = accession, genome = ['cmv', 'spikein', 'human']),
        expand("results/fastqc/{sample}_{direction}_fastqc.html", sample = accession, direction = ["1", "2"]),
        expand("results/fastqc/{sample}_{direction}_processed_fastqc.html", sample = accession, direction = ["1", "2"])
    output:
        "results/multiqc/multiqc_report.html"
    conda:
        "../envs/proseq-qc.yml"
    log:
        out = "log/multiqc.out",
        err = "log/multiqc.err"
    shell:
        "multiqc . --ignore '.snakemake/*' --outdir 'results/multiqc'"
