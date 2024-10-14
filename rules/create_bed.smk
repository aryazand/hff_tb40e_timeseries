rule bam_to_bed:
    input:
        "results/aligned_reads/{sample}_{genome}_extract.bam"
    output:
        "results/bed/{sample}_{genome}.bed"
    conda:
        "../envs/create-bed.yml"
    threads: 5
    log:
        out = "log/bam_to_bed.{sample}_{genome}.out",
        err = "log/bam_to_bed.{sample}_{genome}.err"
    shell:
        "bedtools bamtobed -i {input} > {output} 2> {log.err}"


