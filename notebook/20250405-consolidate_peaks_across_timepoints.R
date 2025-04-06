library(rtracklayer)
library(tidyverse)

###########
# Load bed files with peaks for each time point
###########

files <- list.files("results/peaks", full.names = T)
filenames.split <- files |> 
  strsplit(split = "[/\\.]") |> 
  unlist() |> 
  matrix(ncol = 4, byrow = T)

peaks <- map(files, import.bed)
names(peaks) <- filenames.split[,3]
peaks <- GRangesList(peaks)

#############
# Consolidate peaks across time points
#############

# resize peaks to just the centerpoint
peaks <- resize(peaks, width = 1, fix = "center")
peaks_combined <- unlist(peaks)
peaks_combined$name <- peaks_combined@ranges@NAMES

# consolidate peaks within 30 bp of each other as the same peak, this allows 
# peaks accross timepoints that are slighlty shifted considered the "same" TSS
peaks_combined2 <- peaks_combined |>
  GenomicRanges::reduce(min.gapwidth = 30, with.revmap = T)

revmap <- mcols(peaks_combined2)$revmap
peaks_combined3 <- relist(peaks_combined[unlist(revmap)], revmap) |> 
  map(~GRanges(
    seqnames = seqnames(.x)[1],
    ranges = IRanges(min(start(.x)), max(end(.x))),
    strand = strand(.x)[1],
    mcols = mcols(.x) |> complete()
  ))


