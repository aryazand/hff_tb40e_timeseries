rule trim_reads_pe:
    # Trim adaptors sequences 
    # Use trim_galore to autodetect adaptor sequences 
    input:
        r1="data/fastq/{sample}_1.fastq.gz",
        r2="data/fastq/{sample}_2.fastq.gz"
    output:
        r1=os.path.join(TRIMMED_DIR, "{sample}_1_trimmed.fastq.gz"),
        r2=os.path.join(TRIMMED_DIR, "{sample}_2_trimmed.fastq.gz"),
        report_r1=os.path.join(QC_DIR_TRIMMING, "{sample}_1_trimming_report.txt"),
        report_r2=os.path.join(QC_DIR_TRIMMING, "{sample}_2_trimming_report.txt")
    conda:
        "../envs/trim-galore.yml"
    threads: 5
    log:
        out = "log/trim_reads_{sample}.out",
        err = "log/trim_reads_{sample}.err"
    shell:
	    """
        # make directory if it does not exist
        mkdir -p $(dirname {output.r1})
        mkdir -p $(dirname {output.report_r1})

        # run trim galore
        trim_galore --paired --cores {threads} --gzip \
            --quality 0 \
            --adapter TGGAATTCTCGG --adapter2 GATCGTCGGACT \
            --output_dir $(dirname {output.r1}) \
            {input.r1} {input.r2} \
            2> {log.err} 1> {log.out}
        
        # rename files
        mv $(dirname {output.r1})/{wildcards.sample}_1_val_1.fq.gz {output.r1}
        mv $(dirname {output.r2})/{wildcards.sample}_2_val_2.fq.gz {output.r2}

        # move trimming reports
        mv $(dirname {output.r1})/{wildcards.sample}_1.fastq.gz_trimming_report.txt {output.report_r1}
        mv $(dirname {output.r1})/{wildcards.sample}_2.fastq.gz_trimming_report.txt {output.report_r2}
        """

rule extract_umi: 
    input: 
        r1=os.path.join(TRIMMED_DIR, "{sample}_1_trimmed.fastq.gz"),
        r2=os.path.join(TRIMMED_DIR, "{sample}_2_trimmed.fastq.gz")
    output:
        r1=os.path.join(TRIMMED_DIR, "{sample}_1_trimmed_umi.fastq.gz"),
        r2=os.path.join(TRIMMED_DIR, "{sample}_2_trimmed_umi.fastq.gz"),
        report = os.path.join(UMI_EXTRACT_DIR, "{sample}_umi_report.txt")
    log:
        err = "log/extract_umi_{sample}.err"
    conda: 
        "../envs/umitools.yml"
    shell:
        """
        umi_tools extract -I {input.r1} --bc-pattern=NNNNNNNN --read2-in={input.r2} --stdout={output.r1} --read2-out={output.r2} 2> {log.err} 1> {output.report}
        """
