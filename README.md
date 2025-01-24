# proseq_pipeline
A custom pipeline to analyze proseq data

## How to make apptainer image from container.def file
The container.def file was made as follows, which requires snakemake and spython to generate. It takes all of conda environment files to build a container definition file

```
snakemake --containerize > Dockerfile
spython recipe Dockerfile container.def 
```

To then create an apptainer image, we use the following: 

```
apptainer build container.sif container.def 
```

## UCSC Trackhub
cmv trackhub: <https://genome.ucsc.edu/cgi-bin/hgTracks?genome=KF297339.1&hubUrl=https://raw.githubusercontent.com/aryazand/proseq_pipeline/refs/heads/main/results/UCSCGenomeBrowser/cmv/hub.txt>

human trackhub: <https://genome.ucsc.edu/cgi-bin/hgTracks?genome=GRCh38&hubUrl=https://raw.githubusercontent.com/aryazand/proseq_pipeline/refs/heads/main/results/UCSCGenomeBrowser/human/hub.txt>

## TO DO
1. create recipe to download genomes 
2. create recipe to download fastq 
3. create recipe to get chrom.size data. Ensure chromosome names match UCSC bedgraphtobigwig  