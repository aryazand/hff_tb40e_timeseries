
configfile: "config.yml"

#############################
# Import Snakemake modules
#############################

module proseq_align:
    snakefile: 
        github("aryazand/proseq_align", branch = "main")
    config: config

use rule * from proseq_align exclude sra_to_fastq, gzip_fastq as proseq_align_* 

######################
# Define Rules
#####################

rule all:
    input:
        rules.proseq_align_all.input
    default_target: True