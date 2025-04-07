library(rtracklayer)
library(tidyverse)
library(ggalign)

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
peaks_combined$sample <- factor(peaks_combined@ranges@NAMES, 
                              levels = c("mock", "4hpi", "12hpi", "24hpi", "48hpi", "72hpi"))
peaks_combined$name <- seq(length(peaks_combined))
mcols(peaks_combined) <- mcols(peaks_combined) |> 
  data.frame() |> 
  pivot_wider(id_cols = name, names_from = sample, values_from = score) |> 
  select(-name)

# consolidate peaks within 30 bp of each other as the same peak, this allows 
# peaks accross timepoints that are slighlty shifted considered the "same" TSS
peaks_combined2 <- peaks_combined |>
  GenomicRanges::reduce(min.gapwidth = 30, with.revmap = T)

revmap <- mcols(peaks_combined2)$revmap
mcols(peaks_combined2) <- relist(peaks_combined[unlist(revmap)], revmap) |> 
  map(~mcols(.x) |> data.frame() |> summarise_all(sum, na.rm=T)) |> 
  bind_rows() |> 
  select(mock, 
         `4hpi` = "X4hpi", 
         `12hpi` = "X12hpi",
         `24hpi` = "X24hpi", 
         `48hpi` = "X48hpi",
         `72hpi` = "X72hpi")

txdb <- txdbmaker::makeTxDbFromGFF("data/KF297339.1.gff3", format = "gff3")
cmv.genes <- GenomicFeatures::transcripts(txdb)
nearest(peaks_combined2, cmv.genes)
distanceToNearest(peaks_combined2, cmv.genes)

tss_scores.normalized <- tss_scores |> 
  filter(!(tss %in% c("94775,1,-", "94733,1,+"))) |>
  mutate_at(2:7, .funs = ~ 1000*.x/sum(.x))

tss_scores.top <- tss_scores.normalized |> 
  mutate(sumScore = `4hpi` + `12hpi` + `24hpi` + `48hpi` + `72hpi`) |>
  top_n(100, sumScore)


