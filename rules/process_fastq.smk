rule run_trim:
    # Trim adaptors sequences 
    # Use trim_galore to autodetect adaptor sequences 
    # Don't gzip since downstream application won't work
    input:
        "data/fastq/{id}_1.fastq.gz",
        "data/fastq/{id}_2.fastq.gz"
    output:
        temp("processed_data/fastq/{id}_1_val_1.fq"),
        temp("processed_data/fastq/{id}_2_val_2.fq")
    conda:
        "../envs/process-fastq.yml"
    threads: 5
    log:
        out = "log/run_trim_{id}.out",
        err = "log/run_trim_{id}.err"
    shell:
	    """
        trim_galore --paired --dont_gzip -j {threads} --output_dir processed_data/fastq {input} 2> {log.err} 1> {log.out}
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
        temp("processed_data/fastq/{id}_2_val_2_dedupped.fq.paired.fq")
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

rule rename_processed_files:
    # Rename file
    input:
        "processed_data/fastq/{id}_{direction}_val_{direction}_dedupped.fq.paired.fq"
    output:
        "processed_data/fastq/{id}_{direction}_processed.fastq"
    shell:
        """
        mv {input} {output}
        """