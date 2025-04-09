rule sra_to_fastq:
    # each sra is converted to 3 files based on paired-end reads: *_1.fastq [pair 1], *_2.fastq [pair 2], and *.fastq [reads without pairs]
    output:
        f1 = temp("data/fastq/{sample}_1.fastq"),
        f2 = temp("data/fastq/{sample}_2.fastq")
    conda:
        "../envs/sra-to-fastq.yml"
    threads: 10
    log:
        out = "log/sra_to_fastq_{sample}.out",
        err = "log/sra_to_fastq_{sample}.err"
    params:
        srr=lambda wildcards: sample_table[sample_table.sample_name == wildcards.sample]["srr"].iloc[0]
    threads: 10
    shell:
        """
        # create folder if it does not exist
        mkdir -p data/fastq

        # download SRA data 
        fasterq-dump --split-files -O data/fastq {params.srr} --threads {threads}

        # rename
        mv data/fastq/{params.srr}_1.fastq {output.f1}
        mv data/fastq/{params.srr}_2.fastq {output.f2}
        """

rule gzip_fastq:
    input: 
        f1 = "data/fastq/{sample}_1.fastq",
        f2 = "data/fastq/{sample}_2.fastq"
    output:
        f1 = "data/fastq/{sample}_1.fastq.gz",
        f2 = "data/fastq/{sample}_2.fastq.gz"
    shell:
        """
        gzip {input.f1}
        gzip {input.f2}
        """