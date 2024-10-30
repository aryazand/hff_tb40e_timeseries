rule run_trim:
    # Trim adaptors sequences 
    # Use trim_galore to autodetect adaptor sequences 
    # Don't gzip since downstream application won't work
    input:
        "data/fastq/{id}_1.fastq",
        "data/fastq/{id}_2.fastq"
    output:
        temp("processed_data/fastq/{id}_1_val_1.fq"),
        temp("processed_data/fastq/{id}_2_val_2.fq"),
        report_1 = "results/QC/trimming_reports/{id}_1.fastq_trimming_report.txt",
        report_2 = "results/QC/trimming_reports/{id}_2.fastq_trimming_report.txt"
    conda:
        "../envs/process-fastq.yml"
    threads: 5
    log:
        out = "log/run_trim_{id}.out",
        err = "log/run_trim_{id}.err"
    shell:
	    """
        trim_galore --paired --dont_gzip -j {threads} --output_dir processed_data/fastq {input} 2> {log.err} 1> {log.out}
        mv processed_data/fastq/{wildcards.id}_1.fastq_trimming_report.txt {output.report_1}
        mv processed_data/fastq/{wildcards.id}_2.fastq_trimming_report.txt {output.report_2}
        """

rule run_dedup:
    # Deduplicate reads using seqkit using sequence identity 
    # This is problematic since it doesn't take into account paired-end 
    input:
        "processed_data/fastq/{id}_{direction}_val_{direction}.fq"
    output:
        temp("processed_data/fastq/{id}_{direction}_val_{direction}_dedupped.fq")
    conda:
        "../envs/process-fastq.yml"
    threads: 5
    log:
        out = "log/run_dedup.{id}_{direction}.out",
        err = "log/run_dedup.{id}_{direction}.err"
    shell:
	    """
        seqkit rmdup --threads {threads} -s -o {output} {input} 2> {log.err} 1> {log.out}
        """

rule run_pair:
    # Ensure every read has a pair
    input:
        "processed_data/fastq/{id}_1_val_1_dedupped.fq",
        "processed_data/fastq/{id}_2_val_2_dedupped.fq"
    output:
        temp("processed_data/fastq/{id}_1_val_1_dedupped.fq.paired.fq"),
        temp("processed_data/fastq/{id}_2_val_2_dedupped.fq.paired.fq"),
        temp("processed_data/fastq/{id}_1_val_1_dedupped.fq.single.fq"),
        temp("processed_data/fastq/{id}_2_val_2_dedupped.fq.single.fq")
    conda:
        "../envs/process-fastq.yml"
    threads: 5
    log:
        out = "log/run_pair.{id}.out",
        err = "log/run_pair.{id}.err"
    shell:
        """
        fastq_pair {input} 2> {log.err} 1> {log.out}
        """

rule gzip_processed_fastq:
    # Rename file
    input:
        "processed_data/fastq/{id}_{direction}_val_{direction}_dedupped.fq.paired.fq"
    output:
        "processed_data/fastq/{id}_{direction}_processed.fastq.gz"
    log:
        out = "log/rename_processed_files.{id}_{direction}.out",
        err = "log/rename_processed_files.{id}_{direction}.err"
    shell:
        """
        gzip -cvf {input} > {output}
        rm {input}
        """