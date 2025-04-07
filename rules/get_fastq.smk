rule sra_to_fastq:
    # each sra is converted to 3 files based on paired-end reads: *_1.fastq [pair 1], *_2.fastq [pair 2], and *.fastq [reads without pairs]
    # subsample the fastq_sample option in config file is not empty
    output:
        f1 = temp("data/fastq/{sra}_1.fastq"),
        f2 = temp("data/fastq/{sra}_2.fastq")
    conda:
        "../envs/get-fastq.yml"
    threads: 10
    log:
        out = "log/sra_to_fastq_{sra}.out",
        err = "log/sra_to_fastq_{sra}.err"
    shell:
        """
        fasterq-dump --split-files -O data/fastq {input}
        """
