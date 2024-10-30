rule fastqc:
    input: 
        "data/fastq/{sample}_{direction}.fastq"
    output: 
        "results/QC/fastqc/{sample}_{direction}_fastqc.html"
    conda:
        "../envs/proseq-qc.yml"
    threads: 5
    log:
        out = "log/fastqc_{sample}_{direction}.out",
        err = "log/fastqc_{sample}_{direction}.err"
    shell: 
        "fastqc {input} -t {threads} -o results/QC/fastqc 2> {log.err} 1> {log.out}"


rule fastqc_on_processed_fastq:
    input: 
        "processed_data/fastq/{sample}_{direction}_processed.fastq.gz"
    output: 
        "results/QC/fastqc/{sample}_{direction}_processed_fastqc.html"
    conda:
        "../envs/proseq-qc.yml"
    threads: 5
    log:
        out = "log/fastqc_{sample}_{direction}_processed.out",
        err = "log/fastqc_{sample}_{direction}_processed.err"
    shell: 
        "fastqc {input} -t {threads} -o results/QC/fastqc 2> {log.err} 1> {log.out}"

rule run_multiqc:
    input:
        expand("results/QC/alignment_reports/{sample}.txt", sample = accession),
        expand("results/QC/trimming_reports/{sample}_{direction}.fastq_trimming_report.txt", sample = accession, direction = ["1", "2"]),
        expand("results/QC/fastqc/{sample}_{direction}_fastqc.html", sample = accession, direction = ["1", "2"]),
        expand("results/QC/fastqc/{sample}_{direction}_processed_fastqc.html", sample = accession, direction = ["1", "2"])
    output:
        "results/QC/multiqc/multiqc_report.html"
    conda:
        "../envs/proseq-qc.yml"
    log:
        out = "log/multiqc.out",
        err = "log/multiqc.err"
    shell:
        "multiqc . --ignore '.snakemake/*' --outdir 'results/QC/multiqc' --force"
