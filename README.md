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

## TO DO
1. create recipe to download genomes 
2. create recipe to download fastq 
3. create recipe to get chrom.size data. Ensure chromosome names match UCSC bedgraphtobigwig  