rule align_reads:
    input:
        f1 = "processed_data/fastq/{sample}_1_processed.fastq",
        f2 = "processed_data/fastq/{sample}_2_processed.fastq"
    output:
        temp("results/aligned_reads/{sample}.sam")
    log:
        out = "log/align_reads_{sample}.out",
        err = "log/align_reads_{sample}.err"
    conda:
        "../envs/map-reads.yml"
    threads: 10
    params:
        BOWTIE_INDEX = config['bowtie2']['index'],
        UMI_SIZE = config['bowtie2']['umi_size'],
        PAIRED = config['bowtie2']['paired_options'],
        ADDITIONAL = config['bowtie2']['additional_params']
    shell:
        "bowtie2 -x {params.BOWTIE_INDEX} --threads {threads} --trim5 {params.UMI_SIZE} --trim3 {params.UMI_SIZE} {params.PAIRED} {params.ADDITIONAL} -1 {input.f1} -2 {input.f2} -S {output} 2> {log.err} 1> {log.out}"

rule sam_to_bam:
    input:
        "results/aligned_reads/{sample}.sam"
    output:
        temp("results/aligned_reads/{sample}_allgenomes.bam")
    conda:
        "../envs/map-reads.yml"
    threads: 10
    log:
        out = "log/sam_to_bam.{sample}.out",
        err = "log/sam_to_bam.{sample}.err"
    shell:
        "samtools view --threads {threads} -u {input} | samtools sort --threads 8 -o {output}"

rule sort_bam:
    # Sorts a mapped bam file in preparation for indexing.
    input:
        "results/aligned_reads/{sample}_allgenomes.bam"
    output:
        "results/aligned_reads/{sample}_allgenomes_sorted.bam"
    conda:
        "../envs/map-reads.yml"
    threads: 10
    log:
        out = "log/sort_bam.{sample}_allgenomes.out",
        err = "log/sort_bam.{sample}_allgenomes.err"
    shell:
        """
        samtools sort -@ {threads} {input} -o {output}
        """

rule index_bam:
    input:
        "results/aligned_reads/{sample}_{genome}.bam"
    output:
        "results/aligned_reads/{sample}_{genome}.bam.bai"
    conda:
        "../envs/map-reads.yml"
    threads: 10
    log:
        out = "log/index_bam.{sample}_{genome}.out",
        err = "log/index_bam.{sample}_{genome}.err"
    shell:
       "samtools index {input} 2> {log.err} 1> {log.out}"

rule extract_genome_bam:
    input:
        bam = "results/aligned_reads/{sample}_allgenomes_sorted.bam",
        bai = "results/aligned_reads/{sample}_allgenomes_sorted.bam.bai"
    output:
        "results/aligned_reads/{sample}_{genome}_extract.bam"
    conda:
        "../envs/map-reads.yml"
    threads: 10
    log:
        out = "log/extract_{genome}_bam.{sample}.out",
        err = "log/extract_{genome}_bam.{sample}.err"
    params:   
        genometoextract = lambda wildcards: config['genomes'][wildcards.genome]['pattern_match']
    shell:
        "samtools view --expr 'rname =~ \"{params.genometoextract}\"' --threads {threads} -b -X {input.bam} {input.bai} -o {output} 2> {log.err} 1> {log.out}"

