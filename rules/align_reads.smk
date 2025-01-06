rule create_bowtie2_index:
    input: 
        "data/genome/{combined_species_names}.fna"
    output: 
        multiext("data/genome/{combined_species_names}", 
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    conda:
        "../envs/map-reads.yml"
    log:
        out = "log/create_bowtie2_index_{combined_species_names}.out",
        err = "log/create_bowtie2_index_{combined_species_names}.err"
    threads: 10
    shell:
        """
        bowtie2-build --threads {threads} {input} data/genome/{wildcards.combined_species_names} 2> {log.err} 1> {log.out}
        """

rule align_reads:
    input:
        f1 = "processed_data/fastq/{sample}_R1_processed.fastq.gz",
        f2 = "processed_data/fastq/{sample}_R2_processed.fastq.gz",
        bowtie_index = multiext(f"data/genome/{'_'.join(x for x in SPECIES + SPIKEIN_SPECIES)}", 
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    output:
        sam = temp("results/aligned_reads/{sample}.sam"),
        metrics = "results/QC/alignment_reports/{sample}.txt"
    log:
        out = "log/align_reads_{sample}.out",
        err = "log/align_reads_{sample}.err"
    conda:
        "../envs/map-reads.yml"
    threads: 10
    params:
        BOWTIE2_INDEX = f"data/genome/{'_'.join(x for x in SPECIES + SPIKEIN_SPECIES)}",
        UMI_SIZE = config['bowtie2']['umi_size'],
        PAIRED = config['bowtie2']['paired_options'],
        ADDITIONAL = config['bowtie2']['additional_params']
    shell:
        """
        bowtie2 -x {params.BOWTIE2_INDEX} --threads {threads} --trim5 {params.UMI_SIZE} --trim3 {params.UMI_SIZE} {params.PAIRED} {params.ADDITIONAL} --met-file {output.metrics} -1 {input.f1} -2 {input.f2} -S {output.sam} 2> {log.err} 1> {log.out}
        """

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
        bai = "results/aligned_reads/{sample}_allgenomes_sorted.bam.bai",
        chrom_sizes_file = lambda wc: "data/genome/{species}/" + GENOMES[wc.species] + ".chrom.sizes"
    output:
        "results/aligned_reads/{sample}_{species}_extract.bam"
    conda:
        "../envs/map-reads.yml"
    threads: 10
    log:
        out = "log/extract_{species}_bam.{sample}.out",
        err = "log/extract_{species}_bam.{sample}.err"
    wildcard_constraints:
        genome = '|'.join(GENOMES.keys())
    params:   
        pattern_match = lambda wc, input: '|'.join(x for x in list(pd.read_table(input.chrom_sizes_file, delimiter="\t", header=None).iloc[:,0]))
    shell:
        "samtools view --expr 'rname =~ \"{params.pattern_match}\"' --threads {threads} -b -X {input.bam} {input.bai} -o {output} 2> {log.err} 1> {log.out}"
