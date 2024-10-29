rule bam_to_bedgraph:
    # Create bedgraph file from bam file 
    # Only keep chromosomes from the specific species being processed  
    input:
        bam = "results/aligned_reads/{sample}_{genome}_extract.bam",
        bai = "results/aligned_reads/{sample}_{genome}_extract.bam.bai"
    output:
        bg = "results/tracks/{sample}_{genome}_{direction}.bg"
    log:
        out = "log/bam_to_bedgraph.{sample}_{genome}_{direction}.out",
        err = "log/bam_to_bedgraph.{sample}_{genome}_{direction}.err"
    conda:
        "../envs/create-track.yml"
    threads: 10
    params:
        strand = lambda wildcards: "forward" if wildcards.direction == "for" else "reverse",
        binsize = config['create_track']['binsize'],
        genome_pattern_identifier = lambda wildcards: config['genomes'][wildcards.genome]['pattern_match'],
        track_definition_line = config['create_track']['bedgraph_definition_line']
    shell:
      """
      bamCoverage --outFileFormat bedgraph -p {threads} --binSize {params.binsize} -b {input.bam} --filterRNAstrand {params.strand} -o {output.bg} 2> {log.err} 1> {log.out}
      gawk -i inplace -e '$1 ~ /{params.genome_pattern_identifier}/ {{print $0}}' {output.bg} 2> {log.err} 1> {log.out}
      sed -i '1s/^/{params.track_definition_line}\\n/' {output.bg} 2> {log.err} 1> {log.out}
      """

rule bam_to_bigwig:
    # Create bigwig 
    input:
        bam = "results/aligned_reads/{sample}_{genome}_extract.bam",
        bai = "results/aligned_reads/{sample}_{genome}_extract.bam.bai"
    output:
        bw = "results/tracks/{sample}_{genome}_{direction}.bw"
    log:
        out = "log/bam_to_bigwig.{sample}_{genome}_{direction}.out",
        err = "log/bam_to_bigwig.{sample}_{genome}_{direction}.err"
    threads: 10
    conda:
        "../envs/create-track.yml"
    params:
        binsize = config['create_track']['binsize'],
        strand = lambda wildcards: "forward" if wildcards.direction == "for" else "reverse"
    shell:
      """
      bamCoverage -p {threads} --binSize {params.binsize} -b {input.bam} --filterRNAstrand {params.strand} -o {output.bw} 2> {log.err} 1> {log.out}
      """
