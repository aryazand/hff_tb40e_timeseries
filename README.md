---
bibliography: references.bib
output:
  pdf_document: default
  html_document:
    df_print: paged
---

# Identifying clusters of transcription start regions with similar activation patterns in the cytomegalovirus genome

**Author: Arya Zandvakili**

## Introduction

Cytomegalovirus (CMV) is a beta herpesvirus that infects majority of individuals by age 40 (\~50-80% depending on socioeconomic factors[@colugnati2007]) typically causing mild or no symptoms. After initial infection, CMV remains latent within cells (as an episome) for a person's lifespan. While nearly any cell type can be infected, cells of the myeloid lineage are the primary reservoir for latency. Reactivation of latent CMV in immunosuppressed individuals (AIDS, post-transplant, chemotherapy) can result in severe infections, such as encephalitis, ocular infections, pneumonitis, and colitis. In transplanted individuals, re-activation can result in loss of transplanted organ.

After a cell is infected by CMV, viral genes expression follows a cascade (times related to data from fibroblast cultures)

1.  *Immediate Early* (IE): Genes expressed prior to *de novo* viral protein synthesis (transcribed \< 4 hours post-infection, hpi)

2.  *Early*: Genes that are only expressed after viral protein synthesis initiates, but before viral DNA replication (transcribed \~4-17 hrs hpi)

3.  *Late*: Genes that are only expressed after viral DNA expression (\~48 hrs hpi)

While this pattern of expression is well-established at the level of protein and stable RNA, it is unknown if all transcriptional start regions (TSRs) of the genome have transcriptional activity that can be categorized into “neat” buckets of IE, Early or Late. TSRs encompass all regions of the genome that can initiate transcription (not just promoters of protein-coding gene), such areas of the genome producing non-coding RNA or unstable RNA species (e.g. enhancer). Understanding how transcription is regulated from the CMV genome may not only provide insights to biology of CMV, but may help inform therapeutics.

Here we will use nascent RNA sequencing (Precision Run-On sequencing, PRO-seq) at different time points of CMV infection in cultured human foreskin fibroblasts (HFF) to test the hypothesis that transcription globally from the CMV genome follows the canonical model of IE, Early and Late expression. Additionally, I aim to identify clusters of TSRs that follow similar expression patterns, as this will allow futures studies to determine whether the TSRs in each cluster are regulated by a common mechanism.

## Target Figure

Below is a mock-up of the figure I'd like to produce. It gives an example of what I'd expect if transcriptional activity throughout the genome strongly adhered the canonical model of IE, Early and Late expression. If global transcriptional activity doesn't follow the canonical model we may expect other gene expressions patterns (e.g. areas of the genome producing high levels of activity throughout infection or high levels of expression in 2 time periods).

![](figure_example.png)

## Materials and Methods

### PRO-seq

Briefly, the PRO-seq method exposes cells to biotinylated nucleotides which are then incorporated into nascent RNA strands and resulting in termination of RNA elongation[@mimoso2023]. The biotinylated nucleotides are then isolated with strepavidin pull down to then be sequenced. Treatment of cells with flavopiridol (Flavo) halts Pol II in the promoter-proximal pause zone, resulting in shorter nascent RNA that are less likely to fragment. This allows for more reliable identification of the TSRs (the 5’ end of sequenced RNA is more likely to represent the TSR when the RNA molecules are not fragmented) [@ball2019].

### Production of Raw Data

I will be using the previously published PRO-seq dataset produced by our lab [@ball2022]. This dataset is derived from *in vitro* infection of HFFs with CMV and then PRO-seq being performed at 0, 4, 12-, 24-, 48-, and 72-hours post-infection +/- Flavo at each time point. This is a model of active infection and not latency (as HFF do not harbor latent infection). Prior isolated nascent RNA, a spike-in control of a known number of moth (*Spodoptera*) cells was included with each sample. This allows for comparison expression across time-points at each TSR (i.e. to normalize for difference in read-depth between time points). After isolation of nascent RNA, unique molecular identifier (UMI) sequence adapters are ligated to the ends of the RNA, followed by reverse transcription to cDNA and amplification. cDNA was paired-end sequenced on an Illumina platform. Raw data is FASTQ files representing sequences of nascent RNA produced from CMV, human, and moth genomes.

### Identifying and quantifying TSRs

I will create a snakemake pipeline to do the following:

#### A. Process and Align FASTQ reads

FASTQ will processed as follows: (1) the sequences are deduplicated based on UMI adapters (using `seqkit)`to remove amplification artifacts; (2) adaptors will be trimmed (`trim_galore)`; and (3) any fastq fragment missing its pair will be eliminated (using `fastq_pair`). Reads will be aligned to acombined genome of cytomegalovirus, human, and moth genomes using `bowtie2`.

#### B. Identify TSRs

Create a bed file of just the 5' ends of each mapped read to identify (putative) transcription-start sites (TSSs). TSSs can then be grouped into TSRs. There are many published algorithms to grouping TSSs into TSRs. However, there is no gold-standard. Example algorithms include:

-   `TSRFinder` [@parida2019]: This algorithm analyzes every strand-specific window (size *w*) in the genome and measure the total number of 5' ends mapped within the each window. A window is called a TSR if it (1) meets a user-defined threshold for 5' ends; and (2) does not overlap with a better scoring window. *This is the algorithm my lab has tended to use in the past, but it has weaknesses, such as requiring all TSRs to be of fixed width.*

-   `TSSr` [@lu2021]: This R-package uses the following process to identify TSRs: (1) filter out TSS that don't have enough mapped 5' ends by either user-defined threshold or poisson test compared to background number of 5' ends mapped. (2) Analyze every window (size *w,* default *w* = 100) and within each window, the TSS with the highest 5' ends is identified as the peak. (3) The surrounding weaker TSSs (up to *n* bp away, default *n* = 30) are clustered with the peak to create a TSR. The positions of the 10th to 90th quantiles of TSS signals are defined as the 5' and 3' boundaries of the TSR. The minimum distance between peaks and maximum width of TSRs can be defined.

-   `iTiSS` [@jürges2021]: This algorithm calls TSSs, not TSRs. TSSs are called if they have more 5' ends compared to the average number of 5' ends in surrounding window (size $w$), as implemented as an $F$ score:

    $$
    F = \frac{C_{center} + 1}{\frac{1}{w}\sum_{i=1}^{i=w} C_i + 1} 
    $$

    $C_i$ is the number of 5' ends mapping to position $i$ . Note a pseudocount of 1 is adding to both the numerator and denominator to prevent dividing by zero. A position is called a transcription start site (TSS) if the $F$ score exceeds a certain threshold (which can be user-defined or determined automatically based on the distribution of $F$ scores).

Unfortunately many of these software, including `TSRFinder` and `iTiSS,` are not well-established packages and are not suited for including into snakemake pipeline (e.g. expect software dependencies to be in certain location). I am in the process of creating an R-package (<https://github.com/aryazand/tsrDetectR>) that implements these algorithms. The package will allow me easily determine how much my results change dependent on the algorithm used. It does however have a start-up cost because I'm not experienced in created in R-packages.

#### C. Cluster TSRs based on temporal expression pattern and determine if they are follow the canonical model of CMV gene expression

First, normalize the "signal" of TSSs within the identified TSRs using amount of reads aligned to spike-in moth genome and total number of aligned reads. Second

Second, use hierarchical clustering to cluster TSRs that follow similar temporal gene expression pattern in `R`. Use known IE, Early, and Late genes as comparators to determine if any of the clusters TSRs follow IE, Early, or Late gene expression patterns. *I haven't determined how I will quantify similarity of a TSR (or cluster of TSRs) to a typical IE, Early, or Late gene (or perform statistical testing).*

## Appendix: Reproducing this pipeline

The container.def files describes a container with the exact version of every tool used for this pipeline. You can produce an Apptainer container image from the container.def using the following code:

```         
apptainer build container.sif container.def 
```

Run the pipeline with with the apptainer image using the following command:

```         
snakemake --use-conda --use-apptainer container.sif 
```

## References
