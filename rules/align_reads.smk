rule create_bowtie2_index:
    input: 
        os.path.join(GENOMES_DIR,"{combined_species_names}.fna")
    output: 
        multiext("data/genome/{combined_species_names}", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    conda:
        "../envs/bowtie2.yml"
    log:
        out = "log/create_bowtie2_index_{combined_species_names}.out",
        err = "log/create_bowtie2_index_{combined_species_names}.err"
    threads: 10
    shell:
        """
        bowtie2-build --threads {threads} {input} $(dirname {input}) 2> {log.err} 1> {log.out}
        """

rule align_reads:
    input:
        f1 = os.path.join(TRIMMED_DIR, "{sample}_1_trimmed_umi.fastq.gz"),
        f2 = os.path.join(TRIMMED_DIR, "{sample}_2_trimmed_umi.fastq.gz"),
        bowtie_index = lambda wildcards: multiext(f"data/genome/{"_".join(sample_table[sample_table["sample_name"] == wildcards.sample]["genome_names"].iloc[0].split(", "))}", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    output:
        sam = temp(os.path.join(ALIGNMENT_DIR,"{sample}.sam")),
        metrics = os.path.join(QC_DIR_ALIGNMENT, "{sample}.txt"),
        stats = os.path.join(QC_DIR_ALIGNMENT, "{sample}.stats")
    log:
        out = "log/align_reads_{sample}.out",
    conda:
        "../envs/bowtie2.yml"
    threads: 15
    params:
        BOWTIE_INDEX = lambda wildcards: os.path.join("data/genome/", "_".join(sample_table[sample_table['sample_name'] == wildcards.sample]['genome_names'].iloc[0].split(', '))),
        PAIRED = "--fr --no-discordant",
        ADDITIONAL = ""
    shell:
        """
        bowtie2 -x {params.BOWTIE_INDEX} \
            --threads {threads} \
            {params.PAIRED} {params.ADDITIONAL} \
            --met-file {output.metrics} \
            -1 {input.f1} -2 {input.f2} \
            -S {output.sam} 2> {output.stats} 1> {log.out}
        """

rule sam_to_bam:
    input:
        os.path.join(ALIGNMENT_DIR,"{sample}.sam")
    output:
        temp(os.path.join(ALIGNMENT_DIR, "{sample}_allgenomes.bam"))
    conda:
        "../envs/samtools.yml"
    threads: 10
    log:
        out = "log/sam_to_bam.{sample}.out",
        err = "log/sam_to_bam.{sample}.err"
    shell:
        "samtools view --threads {threads} -u {input} | samtools sort --threads 8 -o {output}"

rule index_bam:
    input:
        os.path.join(ALIGNMENT_DIR, "{sample}_{genome}.bam")
    output:
        os.path.join(ALIGNMENT_DIR, "{sample}_{genome}.bam.bai")
    conda:
        "../envs/samtools.yml"
    threads: 10
    log:
        out = "log/index_bam.{sample}_{genome}.out",
        err = "log/index_bam.{sample}_{genome}.err"
    shell:
       "samtools index {input} 2> {log.err} 1> {log.out}"

rule deduplicate_bam:
    input:
        os.path.join(ALIGNMENT_DIR, "{sample}_allgenomes.bam"),
        os.path.join(ALIGNMENT_DIR, "{sample}_allgenomes.bam.bai")
    output:
        os.path.join(ALIGNMENT_DIR, "{sample}_allgenomes.dedup.bam")
    conda:
        "../envs/umitools.yml"
    threads: 10
    log:
        out = "log/deduplicate_bam.{sample}.out",
        err = "log/deduplicate_bam.{sample}.err"
    params:
        QC_DIR = os.path.join(QC_DIR, "dedup")
    shell:
        """      
        mkdir -p {params.QC_DIR}  
        umi_tools dedup -I {input} --output-stats={params.QC_DIR}/{wildcards.sample} --paired -S {output} 2> {log.err} 1> {log.out}
        """

#rule index_dedup_bam:
#   input:
#       os.path.join(ALIGNMENT_DIR, "{sample}_{genome}.dedup.bam")
#   output:
#       os.path.join(ALIGNMENT_DIR, "{sample}_{genome}.dedup.bam.bai")
#   conda:
#       "../envs/samtools.yml"
#   threads: 10
#   log:
#       out = "log/index_dedup_bam.{sample}_{genome}.out",
#       err = "log/index_dedup_bam.{sample}_{genome}.err"
#   shell:
#      "samtools index {input} 2> {log.err} 1> {log.out}"

rule extract_genome_bam:
    input:
        bam = os.path.join(ALIGNMENT_DIR,"{sample}_allgenomes.dedup.bam"),
        bai = os.path.join(ALIGNMENT_DIR,"{sample}_allgenomes.dedup.bam.bai"),
        chrom_sizes_file = lambda wc: "data/genome/{species}/" + GENOMES[wc.species] + ".chrom.sizes"
    output:
        os.path.join(ALIGNMENT_DIR, "{sample}_{species}.bam")
    conda:
        "../envs/samtools.yml"
    threads: 10
    log:
        out = "log/extract_{species}_bam.{sample}.out",
        err = "log/extract_{species}_bam.{sample}.err"
    wildcard_constraints:
        species = '|'.join(GENOMES.keys())
    params:   
        pattern_match = lambda wc, input: '|'.join(x for x in list(pd.read_table(input.chrom_sizes_file, delimiter="\t", header=None).iloc[:,0]))
    shell:
        "samtools view --expr 'rname =~ \"{params.pattern_match}\"' --threads {threads} -b -X {input.bam} {input.bai} -o {output} 2> {log.err} 1> {log.out}"
