rule fastqc:
    input: 
        os.path.join(RAWFASTQ_DIR, "{sample}_{direction}.fastq.gz")
    output: 
        os.path.join(FASTQC_RAW_DIR, "{sample}_{direction}_fastqc.html")
    wildcard_constraints:
        direction = "[1-2]"
    conda:
        "../envs/proseq-qc.yml"
    threads: 5
    log:
        out = "log/fastqc_{sample}_{direction}.out",
        err = "log/fastqc_{sample}_{direction}.err"
    shell: 
        "fastqc {input} -t {threads} -o $(dirname {output}) 2> {log.err} 1> {log.out}"


rule fastqc_on_processed_fastq:
    input: 
        os.path.join(TRIMMED_DIR, "{sample}_{direction}_val_{direction}_umi.fq.gz")
    output: 
        os.path.join(FASTQC_PROCESSED_DIR, "{sample}_{direction}_val_{direction}_umi_fastqc.html")
    wildcard_constraints:
        direction = "[1-2]"
    conda:
        "../envs/proseq-qc.yml"
    threads: 5
    log:
        out = "log/fastqc_{sample}_{direction}_processed.out",
        err = "log/fastqc_{sample}_{direction}_processed.err"
    shell: 
        "fastqc {input} -t {threads} -o $(dirname {output}) 2> {log.err} 1> {log.out}"

rule run_multiqc:
    input:
        expand(os.path.join(FASTQC_RAW_DIR, "{sample}_{direction}_fastqc.html"), sample = sample_table["sample_name"], direction = ["1","2"]), 
        expand(os.path.join(FASTQC_PROCESSED_DIR, "{sample}_{direction}_val_{direction}_umi_fastqc.html"), sample = sample_table["sample_name"], direction = ["1","2"]), 
    output:
        os.path.join(MULTIQC_DIR, "multiqc_report.html")
    conda:
        "../envs/proseq-qc.yml"
    log:
        out = "log/multiqc.out",
        err = "log/multiqc.err"
    shell:
        "multiqc . --ignore '.snakemake/*' --outdir $(dirname {output}) --force"
