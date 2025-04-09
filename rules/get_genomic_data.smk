rule download_genome:
    # download genome from NCBI
    output: 
        fna = "data/genome/{species}/{accession}.fna",
        gff = "data/genome/{species}/{accession}.gff"
    conda:
        "../envs/get-genome.yml"
    log:
        "log/download_genome_{species}_{accession}.log",
    shell:
        """
        datasets download genome accession {wildcards.accession} --filename {wildcards.species}.zip --include gff3,genome
        unzip {wildcards.species}.zip -d {wildcards.species}
        mv {wildcards.species}/ncbi_dataset/data/{wildcards.accession}/*.fna {output.fna}
        mv {wildcards.species}/ncbi_dataset/data/{wildcards.accession}/*.gff {output.gff}
        sed -i -re 's/(>\\S*)\\s.*/\\1/' {output.fna}
        rm {wildcards.species}.zip
        rm -r {wildcards.species}
        """

rule get_chrom_sizes:
    # get chromosome sizes
    input:
        "data/genome/{species}/{genome}.fna"
    output:
        "data/genome/{species}/{genome}.chrom.sizes"
    conda:
        "../envs/get-genome.yml"
    log:
        "log/get_chrom_sizes_{species}_{genome}.log"
    shell:
        """
        bioawk -cfastx '{{ print $name, length($seq) }}' {input} > {output}
        """

rule concatenate_genomes:
    input:
        expand("data/genome/{species}/{genome}.fna", zip, species = GENOMES.keys(), genome = GENOMES.values())
    output:
        "data/genome/{combined_species_names}.fna"
    log:
        out = "log/concatenate_genomes_{combined_species_names}.out",
        err = "log/concatenate_genomes_{combined_species_names}.err"
    shell:
        """
        cat {input} | sed '/^>/ s/[[:space:]]/\\_/g' > {output}
        """        

rule gff3ToGenePred:
    # convert gff to genePred
    input:
        gff = "data/genome/{species}/{genome}.gff"
    output:
        genePred = "data/genome/{species}/{genome}.genePred"
    conda:
        "../envs/ucsc_tools.yml"
    log:
        "log/gff3ToGenePred_{species}_{genome}.log"
    shell:
        """
        gff3ToGenePred {input.gff} {output.genePred}
        """

rule genePredTobigGenePred:
    # convert genePred to bigGenePred
    input:
        genePred = "data/genome/{species}/{genome}.genePred"
    output:
        biggenePred = "data/genome/{species}/{genome}.biggenePred"
    conda:
        "../envs/ucsc_tools.yml"
    log:
        "log/genePredTobigGenePred_{species}_{genome}.log"
    shell:
        """
        genePredToBigGenePred {input.genePred} {output.biggenePred}
        """

rule create_bigbed:
    # create bigbed for genome features
    input:
        biggenePred = "data/genome/{species}/{genome}.biggenePred",
        chrom_sizes = "data/genome/{species}/{genome}.chrom.sizes"
    output:
        bigbed = "results/UCSCGenomeBrowser/{species}/{genome}/{genome}.bb"
    conda:
        "../envs/ucsc_tools.yml"
    log:
        "log/create_bigbed_{species}_{genome}.log"
    shell:
        """
        sort -k1,1 -k2,2n {input.biggenePred} > {input.biggenePred}.sorted
        wget https://genome.ucsc.edu/goldenPath/help/examples/bigGenePred.as
        bedToBigBed -type=bed12+8 -tab -as=bigGenePred.as {input.biggenePred}.sorted {input.chrom_sizes} {output.bigbed}
        rm bigGenePred.as
        rm {input.biggenePred}.sorted
        """

rule create_2bit:
    # create 2bit file for genome sequence
    input:
        fasta = "data/genome/{species}/{genome}.fna"
    output:
        twobit = "results/UCSCGenomeBrowser/{species}/{genome}/sequence.2bit"
    conda:
        "../envs/ucsc_tools.yml"
    log:
        "log/create_2bit_{species}_{genome}.log"
    shell:
        """
        faToTwoBit {input.fasta} {output.twobit}
        """