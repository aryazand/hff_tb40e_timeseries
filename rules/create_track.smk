rule calculate_readdepth:
    # Calculate normalization factors for tracks
    input:
        bamfiles = expand(["results/aligned_reads/{sample}_{genome}_extract.bam"], sample = accession, genome = ['cmv', 'human', 'spikein']),
        script = "scripts/calculate_normalization_factor.R"
    output:
        "results/QC/total_mapped_reads.tsv"
    log:
        out = "log/calculate_readdepth.out",
        err = "log/calculate_readdepth.err"
    conda:
        "../envs/create-track.yml"
    threads: 10
    shell:
        """
        echo 'Sample\tnum_mapped_pairs' > {output}
        for i in {input.bamfiles}
        do
            sample_name="$(basename $i | cut -d'.' -f1)"
            num_mapped_reads="$(samtools view --threads {threads} --count -f 2 $i)" 
            echo "$sample_name\t$num_mapped_reads" >> {output}
        done

        Rscript {input.script} {output}
        """

rule bam_to_bedgraph:
    # Create bedgraph file from bam file 
    # Only keep chromosomes from the specific species being processed  
    input:
        bam = "results/aligned_reads/{sample}_{genome}_extract.bam",
        bai = "results/aligned_reads/{sample}_{genome}_extract.bam.bai",
        spikein = "results/aligned_reads/{sample}_spikein_extract.bam",
        norm_table = "results/QC/total_mapped_reads.tsv"
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
        sample_name="$(basename {input.bam} | cut -d'_' -f1)"
        scale_factor=$(awk '$1 ~ /'$sample_name'/ {{print $9;}}' {input.norm_table})
        bamCoverage --outFileFormat bedgraph -p {threads} --binSize {params.binsize} --scaleFactor $scale_factor -b {input.bam} --filterRNAstrand {params.strand} -o {output.bg} 2> {log.err} 1> {log.out}
        gawk -i inplace -e '$1 ~ /{params.genome_pattern_identifier}/ {{print $0}}' {output.bg} 2> {log.err} 1> {log.out}
        sed -i '1s/^/{params.track_definition_line}\\n/' {output.bg} 2> {log.err} 1> {log.out}
        """

rule bam_to_bigwig:
    # Create bigwig 
    input:
        bam = "results/aligned_reads/{sample}_{genome}_extract.bam",
        bai = "results/aligned_reads/{sample}_{genome}_extract.bam.bai",
        spikein = "results/aligned_reads/{sample}_spikein_extract.bam",
        norm_table = "results/QC/total_mapped_reads.tsv"
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
        strand = lambda wildcards: "forward" if wildcards.direction == "for" else "reverse",
    shell:
        """
        sample_name="$(basename {input.bam} | cut -d'_' -f1)"
        scale_factor=$(awk '$1 ~ /'$sample_name'/ {{print $9;}}' {input.norm_table})
        bamCoverage -p {threads} --binSize {params.binsize} --scaleFactor $scale_factor -b {input.bam} --filterRNAstrand {params.strand} -o {output.bw} 2> {log.err} 1> {log.out}
        """
