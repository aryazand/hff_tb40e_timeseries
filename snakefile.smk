fastq_files, = glob_wildcards("data/fastq/{fname}_1.fq")

rule all:
    input:
        expand("results/fastqc/{sample}_1_fastqc.html", sample = fastq_files),
        expand("results/fastqc/{sample}_2_fastqc.html", sample = fastq_files),
        expand("results/fastqc/{sample}_1_val_1_dedupped.fq.paired_fastqc.html", sample = fastq_files),
        expand("results/fastqc/{sample}_2_val_2_dedupped.fq.paired_fastqc.html", sample = fastq_files),
        expand("results/aligned_reads/{sample}.sam", sample = fastq_files),

rule fastqc:
    input: 
        "data/fastq/{sample}.fq"
    output: 
        "results/fastqc/{sample}_fastqc.html"
    log:
        out = "log/fastqc_{sample}.out",
        err = "log/fastqc_{sample}.err"
    shell: 
        "fastqc {input} -t 10 -o results/fastqc 2> {log.err} 1> {log.out}"

##########################
# Process fastq files
##########################

rule run_trim:
    input:
        "data/fastq/{sample}_1.fq",
        "data/fastq/{sample}_2.fq"
    output:
        "processed_data/trimmed_fastq/{sample}_1_val_1.fq",
        "processed_data/trimmed_fastq/{sample}_2_val_2.fq"
    log:
        out = "log/run_trim_{sample}.out",
        err = "log/run_trim_{sample}.err"
    shell:
	    "trim_galore --paired --dont_gzip -j 5 --output_dir processed_data/trimmed_fastq {input} 2> {log.err} 1> {log.out}"

rule run_dedup:
    input:
        "processed_data/trimmed_fastq/{sample}_{direction}_val_{direction}.fq"
    output:
        "processed_data/trimmed_dedup_fastq/{sample}_{direction}_val_{direction}_dedupped.fq"
    log:
        out = "log/run_dedup.{sample}_{direction}.out",
        err = "log/run_dedup.{sample}_{direction}.err"
    shell:
	    "seqkit rmdup -s -o {output} {input} 2> {log.err} 1> {log.out}"

rule run_pair:
    input:
        f1 = "processed_data/trimmed_dedup_fastq/{sample}_1_val_1_dedupped.fq",
        f2 = "processed_data/trimmed_dedup_fastq/{sample}_2_val_2_dedupped.fq"
    output:
        f1 = "processed_data/trimmed_dedup_fastq_paired/{sample}_1_val_1_dedupped.fq.paired.fq",
        f2 = "processed_data/trimmed_dedup_fastq_paired/{sample}_2_val_2_dedupped.fq.paired.fq"
    log:
        out = "log/run_pair.{sample}.out",
        err = "log/run_pair.{sample}.err"
    shell:
        """
        fastq_pair {input} 2> {log.err} 1> {log.out}
        mv processed_data/trimmed_dedup_fastq/{wildcards.sample}_1_val_1_dedupped.fq.paired.fq {output.f1}
        mv processed_data/trimmed_dedup_fastq/{wildcards.sample}_2_val_2_dedupped.fq.paired.fq {output.f2}
        """

rule run_fastqc_on_processed_reads:
    input: 
        "processed_data/trimmed_dedup_fastq_paired/{sample}_{direction}_val_{direction}_dedupped.fq.paired.fq"
    output: 
        "results/fastqc/{sample}_{direction}_val_{direction}_dedupped.fq.paired_fastqc.html"
    log:
        out = "log/fastqc_{sample}_{direction}_val_{direction}_dedupped.out",
        err = "log/fastqc_{sample}_{direction}_val_{direction}_dedupped.err"
    shell: 
        "fastqc {input} -t 10 -o results/fastqc 2> {log.err} 1> {log.out}"

##########################
# Align reads
##########################

rule align_reads:
    input:
        f1 = "processed_data/trimmed_dedup_fastq_paired/{sample}_1_val_1_dedupped.fq.paired.fq",
        f2 = "processed_data/trimmed_dedup_fastq_paired/{sample}_2_val_2_dedupped.fq.paired.fq"
    output:
        "results/aligned_reads/{sample}.sam"
    log:
        out = "log/align_reads_{sample}.out",
        err = "log/align_reads_{sample}.err"
    params:
        BOWTIE_INDEX = "data/genome/human_tb40E_spodoptera_combined_genome",
        UMI_SIZE = 8
    shell:
        "bowtie -x {params.BOWTIE_INDEX} --threads 4 --trim5 {params.UMI_SIZE} --trim3 {params.UMI_SIZE} --chunkmbs 500 --fr --best --sam --allow-contain --fullref -1 {input.f1} -2 {input.f2} {output} 2> {log.err} 1> {log.out}"