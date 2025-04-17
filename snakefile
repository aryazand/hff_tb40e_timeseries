
configfile: "config.yml"

#############################
# Import Snakemake modules
#############################

module proseq_align:
    snakefile: 
        github("aryazand/align_proseq", path = "snakefile", branch = "main")
    config: config

use rule * from proseq_align as proseq_align_* 

######################
# Define Rules
#####################

rule all:
    input:
        rules.proseq_align_all.input
    default_target: True