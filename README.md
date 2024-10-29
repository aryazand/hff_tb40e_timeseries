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

