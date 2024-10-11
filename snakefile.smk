fastq_files, = glob_wildcards("data/fastq/{fname}_1.fastq.gz")

rule all:
    input:
        expand("results/fastqc/{sample}_{direction}_fastqc.html", sample = fastq_files, direction = ["1", "2"]),
        expand("results/fastqc/{sample}_{direction}_val_{direction}_dedupped.fq.paired_fastqc.html", sample = fastq_files, direction = ["1", "2"]),
        expand("processed_data/trimmed_dedup_fastq_paired/{sample}_{direction}_val_{direction}_dedupped.fq.paired.fq", sample = fastq_files, direction = ["1", "2"]),
        expand("results/aligned_reads/{sample}_cmv.bed", sample = fastq_files),
        expand("results/aligned_reads/{sample}_{genome}_{direction}.bg", sample = fastq_files, genome = "cmv", direction = ["for", "rev"]),
        expand("results/aligned_reads/{sample}_{genome}.bam", sample = fastq_files, genome = ["cmv", "allgenomes", "spikein"]),
        expand("results/aligned_reads/{sample}_{genome}.bam.bai", sample = fastq_files, genome = ["cmv", "allgenomes", "spikein"])


rule fastqc:
    input: 
        "data/fastq/{sample}.fastq.gz"
    output: 
        "results/fastqc/{sample}_fastqc.html"
    log:
        out = "log/fastqc_{sample}.out",
        err = "log/fastqc_{sample}.err"
    shell: 
        "fastqc {input} -t 10 -o results/fastqc 2> {log.err} 1> {log.out}"

##########################
# Process fastq files
##########################

rule run_trim:
    input:
        "data/fastq/{sample}_1.fastq.gz",
        "data/fastq/{sample}_2.fastq.gz"
    output:
        "processed_data/trimmed_fastq/{sample}_1_val_1.fq",
        "processed_data/trimmed_fastq/{sample}_2_val_2.fq"
    log:
        out = "log/run_trim_{sample}.out",
        err = "log/run_trim_{sample}.err"
    shell:
	    "trim_galore --paired --dont_gzip -j 5 --output_dir processed_data/trimmed_fastq {input} 2> {log.err} 1> {log.out}"

rule run_dedup:
    input:
        "processed_data/trimmed_fastq/{sample}_{direction}_val_{direction}.fq"
    output:
        "processed_data/trimmed_dedup_fastq/{sample}_{direction}_val_{direction}_dedupped.fq"
    log:
        out = "log/run_dedup.{sample}_{direction}.out",
        err = "log/run_dedup.{sample}_{direction}.err"
    shell:
	    """
        seqkit rmdup -s -o {output} {input} 2> {log.err} 1> {log.out}
        rm {input}
        """

rule run_pair:
    input:
        f1 = "processed_data/trimmed_dedup_fastq/{sample}_1_val_1_dedupped.fq",
        f2 = "processed_data/trimmed_dedup_fastq/{sample}_2_val_2_dedupped.fq"
    output:
        f1 = "processed_data/trimmed_dedup_fastq_paired/{sample}_1_val_1_dedupped.fq.paired.fq",
        f2 = "processed_data/trimmed_dedup_fastq_paired/{sample}_2_val_2_dedupped.fq.paired.fq"
    log:
        out = "log/run_pair.{sample}.out",
        err = "log/run_pair.{sample}.err"
    shell:
        """
        fastq_pair {input} 2> {log.err} 1> {log.out}
        mv processed_data/trimmed_dedup_fastq/{wildcards.sample}_1_val_1_dedupped.fq.paired.fq {output.f1}
        mv processed_data/trimmed_dedup_fastq/{wildcards.sample}_2_val_2_dedupped.fq.paired.fq {output.f2}
        rm {input.f1} {input.f2} 
        """

rule run_fastqc_on_processed_reads:
    input: 
        "processed_data/trimmed_dedup_fastq_paired/{sample}_{direction}_val_{direction}_dedupped.fq.paired.fq"
    output: 
        "results/fastqc/{sample}_{direction}_val_{direction}_dedupped.fq.paired_fastqc.html"
    log:
        out = "log/fastqc_{sample}_{direction}_val_{direction}_dedupped.out",
        err = "log/fastqc_{sample}_{direction}_val_{direction}_dedupped.err"
    shell: 
        "fastqc {input} -t 10 -o results/fastqc 2> {log.err} 1> {log.out}"

##########################
# Align reads
##########################

rule align_reads:
    input:
        f1 = "processed_data/trimmed_dedup_fastq_paired/{sample}_1_val_1_dedupped.fq.paired.fq",
        f2 = "processed_data/trimmed_dedup_fastq_paired/{sample}_2_val_2_dedupped.fq.paired.fq"
    output:
        "results/aligned_reads/{sample}.sam"
    log:
        out = "log/align_reads_{sample}.out",
        err = "log/align_reads_{sample}.err"
    params:
        BOWTIE_INDEX = "data/genome/human_tb40E_spodoptera_combined_genome",
        UMI_SIZE = 8
    shell:
        "bowtie -x {params.BOWTIE_INDEX} --threads 10 --trim5 {params.UMI_SIZE} --trim3 {params.UMI_SIZE} --chunkmbs 500 --fr --best --sam --allow-contain --fullref -1 {input.f1} -2 {input.f2} {output} 2> {log.err} 1> {log.out}"

rule sam_to_bam:
    input:
        "results/aligned_reads/{sample}.sam"
    output:
        "results/aligned_reads/{sample}_allgenomes.bam"
    log:
        out = "log/sam_to_bam.{sample}.out",
        err = "log/sam_to_bam.{sample}.err"
    shell:
        "samtools view --threads 8 -u {input} | samtools sort --threads 8 -o {output}"

rule index_bam:
    input:
        "results/aligned_reads/{sample}_{genome}.bam"
    output:
        "results/aligned_reads/{sample}_{genome}.bam.bai"
    log:
        out = "log/index_bam.{sample}_{genome}.out",
        err = "log/index_bam.{sample}_{genome}.err"
    shell:
       "samtools index {input} 2> {log.err} 1> {log.out}"

rule extract_cmv_bam:
    input:
        bam = "results/aligned_reads/{sample}_allgenomes.bam",
        bai = "results/aligned_reads/{sample}_allgenomes.bam.bai"
    output:
        "results/aligned_reads/{sample}_cmv.bam"
    log:
        out = "log/extract_cmv_bam.{sample}.out",
        err = "log/extract_cmv_bam.{sample}.err"
    params:   
        cmv_genome = "KF297339.1"
    shell:
        "samtools view --threads 8 -b -X {input.bam} {input.bai} {params.cmv_genome} -o {output} 2> {log.err} 1> {log.out}"

rule extract_spikein_bam:
    input:
        bam = "results/aligned_reads/{sample}_allgenomes.bam",
        bai = "results/aligned_reads/{sample}_allgenomes.bam.bai"
    output:
        "results/aligned_reads/{sample}_spikein.bam"
    log:
        out = "log/extract_spikein_bam.{sample}.out",
        err = "log/extract_spikein_bam.{sample}.err"
    params:   
        spikein_genome = "Spodoptera"
    shell:
        "samtools view --threads 8 -b --expr 'rname =~ \"{params.spikein_genome}\"' -X {input.bam} {input.bai} -o {output} 2> {log.err} 1> {log.out}"


rule bam_to_bed:
    input:
        "results/aligned_reads/{sample}_cmv.bam"
    output:
        "results/aligned_reads/{sample}_cmv.bed"
    log:
        out = "log/bam_to_bed.{sample}.out",
        err = "log/bam_to_bed.{sample}.err"
    shell:
        "bedtools bamtobed -i {input} > {output} 2> {log.err}"


rule bam_to_bedgraph:
    input:
        bam = "results/aligned_reads/{sample}_{genome}.bam"
    output:
        for_bg = "results/aligned_reads/{sample}_{genome}_for.bg",
        rev_bg = "results/aligned_reads/{sample}_{genome}_rev.bg"
    log:
        out_for = "log/bam_to_bedgraph.{sample}_{genome}_for.out",
        out_rev = "log/bam_to_bedgraph.{sample}_{genome}_rev.out",
        err_for = "log/bam_to_bedgraph.{sample}_{genome}_for.err",
        err_rev = "log/bam_to_bedgraph.{sample}_{genome}_rev.err"
    params:
        threads = 10,
        binsize = 1,
        genome_pattern_identifier = "KF297339",
        track_definition_line = "track type=bedGraph"
    shell:
      """
      bamCoverage --outFileFormat bedgraph -p {params.threads} --binSize {params.binsize} -b {input.bam} --filterRNAstrand forward -o {output.for_bg} 2> {log.err_for} 1> {log.out_for}
      awk -i inplace -e '$1 ~ /{params.genome_pattern_identifier}/ {{print $0}}' {output.for_bg} 2> {log.err_for} 1> {log.out_for}
      sed -i '1s/^/{params.track_definition_line}\\n/' {output.for_bg} 2> {log.err_for} 1> {log.out_for}
      
      bamCoverage --outFileFormat bedgraph -p {params.threads} --binSize {params.binsize} -b {input.bam} --filterRNAstrand reverse -o {output.rev_bg} 2> {log.err_rev} 1> {log.out_rev}
      awk -i inplace -e '$1 ~ /{params.genome_pattern_identifier}/ {{print $0}}' {output.rev_bg} 2> {log.err_rev} 1> {log.out_rev}
      sed -i '1s/^/{params.track_definition_line}\\n/' {output.rev_bg} 2> {log.err_rev} 1> {log.out_rev}
      """