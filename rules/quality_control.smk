rule fastqc:
    input: 
        "processed_data/fastq/{id}_R{direction}" + sample_ext
    output: 
        "results/QC/fastqc/{id}_R{direction}_fastqc.html"
    wildcard_constraints:
        direction = "[1-2]"
    conda:
        "../envs/proseq-qc.yml"
    threads: 5
    log:
        out = "log/fastqc_{id}_{direction}.out",
        err = "log/fastqc_{id}_{direction}.err"
    shell: 
        "fastqc {input} -t {threads} -o results/QC/fastqc 2> {log.err} 1> {log.out}"


rule fastqc_on_processed_fastq:
    input: 
        "processed_data/fastq/{sample}_R{direction}_processed.fastq.gz"
    output: 
        "results/QC/fastqc/{sample}_R{direction}_processed_fastqc.html"
    wildcard_constraints:
        direction = "[1-2]"
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
        expand("results/QC/alignment_reports/{sample}.txt", sample = sample_names),
        expand("results/QC/trimming_reports/{sample}_{direction}{ext}_trimming_report.txt", sample = sample_names, direction = ["R1", "R2"], ext = sample_ext),
        expand("results/QC/fastqc/{sample}_{direction}_fastqc.html", sample = sample_names, direction = ["R1", "R2"]),
        expand("results/QC/fastqc/{sample}_{direction}_processed_fastqc.html", sample = sample_names, direction = ["R1", "R2"])
    output:
        "results/QC/multiqc/multiqc_report.html"
    conda:
        "../envs/proseq-qc.yml"
    log:
        out = "log/multiqc.out",
        err = "log/multiqc.err"
    shell:
        "multiqc . --ignore '.snakemake/*' --outdir 'results/QC/multiqc' --force"


