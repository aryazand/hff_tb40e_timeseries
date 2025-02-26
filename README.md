# Identifying clusters of transcription start regions in cytomegalovirus genome with similiar activation patterns 

Here we will analyze nascent RNA via PRO-seq at different timepoints of cytomegalovirus (CMV) infection of fibroblasts to identifying transcription start regions (TSRs) that follow similar activation patterns. We will analyze the PRO-seq dataset is hosted on NCBI GEO GSE185763, which is associated with publication doi: 10.3390/v14040779.

## Introduction
TODO: Add details

## Description of data
Human foreskin fibroblast were infected with cytomegalovirus and PRO-seq was performed with and without flavopriridol at 0, 4, 12, 24, 48, and 72 hours post-infection. 

TODO: Add more description 


## Reproducing this pipeline

### How to make apptainer image from container.def file
The container.def file was made as follows, which requires snakemake and spython to generate. It takes all of conda environment files to build a container definition file

```
snakemake --containerize > Dockerfile
spython recipe Dockerfile container.def 
```

To then create an apptainer image, we use the following: 

```
apptainer build container.sif container.def 
```


