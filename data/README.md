# README: Metadata regarding raw data


## FASTQ files

The files in the fastq folder are downloaded from NCBI SRA database. They represent PRO-seq data derived from human foreskin fibroblasts (HFF) at different hours post infection (hpi) with cytomegalovirus. The data is associated with NCBI GEO dataset GSE185762. If these fastq files are not present, the snakemake pipeline automatically dowloads SRA usings `prefetch` function from the `sra-toolkit` and converts the SRA to FASTQ using `fasterq-dump`. See recipe rules/get_fastq.smk for exact code to produce this data.

| File                       | Description          | NCBI SRA Accession |
|----------------------------|----------------------|--------------------|
| fastq/SRR16301851.fastq.gz | Mock DMSO PRO-Seq    | SRR16301851        |
| fastq/SRR16301852.fastq.gz | Mock Flavo PRO-Seq   | SRR16301852        |
| fastq/SRR16301853.fastq.gz | 4 hpi DMSO PRO-Seq   | SRR16301853        |
| fastq/SRR16301854.fastq.gz | 4 hpi Flavo PRO-Seq  | SRR16301854        |
| fastq/SRR16301855.fastq.gz | 12 hpi DMSO PRO-Seq  | SRR16301855        |
| fastq/SRR16301856.fastq.gz | 12 hpi Flavo PRO-Seq | SRR16301856        |
| fastq/SRR16301857.fastq.gz | 24 hpi DMSO PRO-Seq  | SRR16301857        |
| fastq/SRR16301858.fastq.gz | 24 hpi Flavo PRO-Seq | SRR16301858        |
| fastq/SRR16301859.fastq.gz | 48 hpi DMSO PRO-Seq  | SRR16301859        |
| fastq/SRR16301860.fastq.gz | 48 hpi Flavo PRO-Seq | SRR16301860        |
| fastq/SRR16301861.fastq.gz | 72 hpi DMSO PRO-Seq  | SRR16301861        |
| fastq/SRR16301862.fastq.gz | 72 hpi Flavo PRO-Seq | SRR16301862        |

## Genomes

Using NCBI `datasets` toolkit, the following genomes are downloaded NCBI. See snakemake recipe rules/get_genomic_data.smk for code. These zip file contain both fasta sequences for the genomes as well as GFF files describing genomic features. 

| File                  | Short Label                         | Genbank Accession | Description                                |
|-----------------------|-------------------------------------|-------------------|--------------------------------------------|
| genome/GRCh38.zip     | human hg38 genome                   | GCF_000001405.40  | human genome hg38                          |
| genome/KF297339.1.zip | human cytomegalovirus TB40-E genome | GCA_027926625.1   | cytomegalovirus genome TB-40E (KF297339.1) |
| genome/spodoptera.zip | Spodoptera genome                   | GCF_023101765.2   | used as spike-in control in experiments    |
