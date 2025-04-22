
configfile: "config.yml"

#############################
# Import Snakemake modules
#############################

module proseq_align:
    snakefile: 
        github("aryazand/align_proseq", path = "snakefile", branch = "main")
    config: config

use rule * from proseq_align as proseq_align_* 

module bam_to_ucschub:
    snakefile: "../bam_to_ucschub/snakefile"
    config: config

use rule * from bam_to_ucschub as bam_to_ucschub_* 

######################
# Define Rules
#####################

rule all:
    input:
        rules.proseq_align_all.input,
        rules.bam_to_ucschub_all.input
    default_target: True

rule multiqc:
    input:
        rules.proseq_align_multiqc.input
