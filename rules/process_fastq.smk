rule trim_reads_pe:
    # Trim adaptors sequences 
    # Use trim_galore to autodetect adaptor sequences 
    input:
        r1="data/fastq/{sample}_1.fastq.gz",
        r2="data/fastq/{sample}_2.fastq.gz"
    output:
        r1=temp(os.path.join(TRIMMED_DIR, "{sample}_1_val_1.fq.gz")),
        r2=temp(os.path.join(TRIMMED_DIR, "{sample}_2_val_2.fq.gz")),
        r1_report=os.path.join(QC_DIR_TRIMMING, "{sample}_1.fastq.gz_trimming_report.txt"),
        r2_report=os.path.join(QC_DIR_TRIMMING, "{sample}_2.fastq.gz_trimming_report.txt")
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
        mkdir -p $(dirname {output.r1_report})

        # run trim galore
        trim_galore --paired -j {threads} --output_dir $(dirname {output.r1}) {input.r1} {input.r2} 2> {log.err} 1> {log.out}
        
        # move trimming reports
        mv $(dirname {output.r1})/{wildcards.sample}_1.fastq.gz_trimming_report.txt {output.r1_report}
        mv $(dirname {output.r2})/{wildcards.sample}_2.fastq.gz_trimming_report.txt {output.r2_report}
        """

rule extract_umi: 
    input: 
        r1=os.path.join(TRIMMED_DIR, "{sample}_1_val_1.fq.gz"),
        r2=os.path.join(TRIMMED_DIR, "{sample}_2_val_2.fq.gz")
    output:
        r1=os.path.join(TRIMMED_DIR, "{sample}_1_val_1_umi.fq.gz"),
        r2=os.path.join(TRIMMED_DIR, "{sample}_2_val_2_umi.fq.gz")
    params:
        QC_DIR = QC_DIR_UMIEXTRACT
    log:
        err = "log/extract_umi_{sample}.err"
    conda: 
        "../envs/umitools.yml"
    shell:
        """
        umi_tools extract -I {input.r1} --bc-pattern=NNNNNNNN --read2-in={input.r2} --stdout={output.r1} --read2-out={output.r2} 2> {log.err} 1> {params.QC_DIR}
        """
