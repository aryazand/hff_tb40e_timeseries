rule concatenate_fastq:
    # Concatenate files from the same sample into one fastq 
    # Files may have been separated into lanes
    input:
       lambda wc: config['fastq']['directory'] + samples[(samples.sample_name == wc.id) & (samples.direction == wc.direction)].filename
    output: 
        temp("processed_data/fastq/{id}_{direction}.fastq.gz")
    wildcard_constraints:
        direction = "R[1-2]"
    log:
        out = "log/concatenate_fastq_{id}_{direction}.out",
        err = "log/concatenate_fastq_{id}_{direction}.err"
    shell:
        """
        echo 'input is {input}'
        cat {input} > {output}
        """

rule concatenate_fastq_uncompressed:
    # Concatenate files from the same sample into one fastq
    # Files may have been separated into lanes
    input:
       lambda wc: config['fastq']['directory'] + samples[(samples.sample_name == wc.id) & (samples.direction == wc.direction)].filename
    output: 
        temp("processed_data/fastq/{id}_{direction}.fastq")
    wildcard_constraints:
        direction = "R[1-2]"
    log:
        out = "log/concatenate_fastq_{id}_{direction}.out",
        err = "log/concatenate_fastq_{id}_{direction}.err"
    shell:
        """
        cat {input} > {output}
        """

rule trim_reads_pe:
    # Trim adaptors sequences 
    # Use trim_galore to autodetect adaptor sequences 
    # Don't gzip since downstream application won't work
    # ext = set([''.join(pathlib.Path(x).suffixes) for x in samples[samples.sample_name == wc.id].filename.to_list()])
    input:
        "processed_data/fastq/{id}_R1" + sample_ext,
        "processed_data/fastq/{id}_R2" + sample_ext
    output:
        temp("processed_data/fastq/{id}_R1_val_1.fq"),
        temp("processed_data/fastq/{id}_R2_val_2.fq"),
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

rule move_trimming_report:
    # Move trimming report ot QC area
    input:
        report_1 = "processed_data/fastq/{id}_R1.{ext}_trimming_report.txt",
        report_2 = "processed_data/fastq/{id}_R2.{ext}_trimming_report.txt"
    output:
        report_1 = "results/QC/trimming_reports/{id}_R1.{ext}_trimming_report.txt",
        report_2 = "results/QC/trimming_reports/{id}_R2.{ext}_trimming_report.txt"
    shell:
        """
        mv {input.report_1} {output.report_1}
        mv {input.report_2} {output.report_2}
        """

rule run_dedup:
    # Deduplicate reads using seqkit using sequence identity 
    # This is problematic since it doesn't take into account paired-end 
    input:
        "processed_data/fastq/{id}_R{direction}_val_{direction}.fq"
    output:
        temp("processed_data/fastq/{id}_R{direction}_val_{direction}_dedupped.fq")
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
        "processed_data/fastq/{id}_R1_val_1_dedupped.fq",
        "processed_data/fastq/{id}_R2_val_2_dedupped.fq"
    output:
        temp("processed_data/fastq/{id}_R1_val_1_dedupped.fq.paired.fq"),
        temp("processed_data/fastq/{id}_R2_val_2_dedupped.fq.paired.fq"),
        temp("processed_data/fastq/{id}_R1_val_1_dedupped.fq.single.fq"),
        temp("processed_data/fastq/{id}_R2_val_2_dedupped.fq.single.fq")
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
        "processed_data/fastq/{id}_R{direction}_val_{direction}_dedupped.fq.paired.fq"
    output:
        "processed_data/fastq/{id}_R{direction}_processed.fastq.gz"
    log:
        out = "log/rename_processed_files.{id}_{direction}.out",
        err = "log/rename_processed_files.{id}_{direction}.err"
    shell:
        """
        gzip -cvf {input} > {output}
        rm {input}
        """