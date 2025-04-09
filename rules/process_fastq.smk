rule trim_reads_pe:
    # Trim adaptors sequences 
    # Use trim_galore to autodetect adaptor sequences 
    input:
        r1="data/fastq/{sample}_1.fastq.gz",
        r2="data/fastq/{sample}_2.fastq.gz"
    output:
        r1=os.path.join(TRIMMED_DIR, "{sample}_1_val_1.fq.gz"),
        r2=os.path.join(TRIMMED_DIR, "{sample}_2_val_2.fq.gz"),
        r1_report=os.path.join(QC_DIR, "{sample}_1.fastq.gz_trimming_report.txt"),
        r2_report=os.path.join(QC_DIR, "{sample}_2.fastq.gz_trimming_report.txt")
    params:
        TRIMMED_DIR=TRIMMED_DIR,
        QC_DIR=os.path.join(QC_DIR, "trimming_reports")
    conda:
        "../envs/trim-galore.yml"
    threads: 5
    log:
        out = "log/trim_reads_{sample}.out",
        err = "log/trim_reads_{sample}.err"
    shell:
	    """
        # make directory if it does not exist
        mkdir -p {params.TRIMMED_DIR}
        mkdir -p {params.QC_DIR}

        # run trim galore
        trim_galore --paired -j {threads} --output_dir {TRIMMED_DIR} {input.r1} {input.r2} 2> {log.err} 1> {log.out}
        
        # move trimming reports
        mv {params.TRIMMED_DIR}/{wildcards.sample}_1.fastq.gz_trimming_report.txt {output.r1_report}
        mv {params.TRIMMED_DIR}/{wildcards.sample}_2.fastq.gz_trimming_report.txt {output.r2_report}
        """