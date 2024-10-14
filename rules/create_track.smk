rule bam_to_bedgraph:
    # Create bedgraph file from bam file 
    # Only keep chromosomes from the specific species being processed  
    input:
        bam = "results/aligned_reads/{sample}_{genome}_extract.bam",
        bai = "results/aligned_reads/{sample}_{genome}_extract.bam.bai"
    output:
        for_bg = "results/tracks/{sample}_{genome}_for.bg",
        rev_bg = "results/tracks/{sample}_{genome}_rev.bg"
    log:
        out_for = "log/bam_to_bedgraph.{sample}_{genome}_for.out",
        out_rev = "log/bam_to_bedgraph.{sample}_{genome}_rev.out",
        err_for = "log/bam_to_bedgraph.{sample}_{genome}_for.err",
        err_rev = "log/bam_to_bedgraph.{sample}_{genome}_rev.err"
    conda:
        "../envs/create-track.yml"
    threads: 10
    params:
        binsize = config['create_bedgraph']['binsize'],
        genome_pattern_identifier = lambda wildcards: config['genomes'][wildcards.genome]['pattern_match'],
        track_definition_line = config['create_bedgraph']['track_definition_line']
    shell:
      """
      bamCoverage --outFileFormat bedgraph -p {threads} --binSize {params.binsize} -b {input.bam} --filterRNAstrand forward -o {output.for_bg} 2> {log.err_for} 1> {log.out_for}
      awk -i inplace -e '$1 ~ /{params.genome_pattern_identifier}/ {{print $0}}' {output.for_bg} 2> {log.err_for} 1> {log.out_for}
      sed -i '1s/^/{params.track_definition_line}\\n/' {output.for_bg} 2> {log.err_for} 1> {log.out_for}
      
      bamCoverage --outFileFormat bedgraph -p {threads} --binSize {params.binsize} -b {input.bam} --filterRNAstrand reverse -o {output.rev_bg} 2> {log.err_rev} 1> {log.out_rev}
      awk -i inplace -e '$1 ~ /{params.genome_pattern_identifier}/ {{print $0}}' {output.rev_bg} 2> {log.err_rev} 1> {log.out_rev}
      sed -i '1s/^/{params.track_definition_line}\\n/' {output.rev_bg} 2> {log.err_rev} 1> {log.out_rev}
      """