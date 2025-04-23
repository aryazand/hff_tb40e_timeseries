library(Rsamtools)
library(tidyverse)
library(GenomicAlignments)

extract_5prime_ends <- function(bam_file, output_file = NULL, strand_specific = TRUE) {
  # Print status
  cat("Processing", bam_file, "...\n")
  
  # Open BAM file
  bamfile <- BamFile(bam_file)
  
  # Read BAM file
  param <- ScanBamParam(
    flag = scanBamFlag(isProperPair = TRUE, isPaired = TRUE),
    what = c("rname", "pos", "qwidth", "strand", "mpos", "isize"))
  reads <- readGAlignments(bamfile, param = param)
  
  # Convert to GRanges
  gr <- granges(reads)

  five_prime_ends <- GRanges(
      seqnames = seqnames(gr),
      ranges = IRanges(
        start = ifelse(strand(gr) == "+", start(gr), end(gr)),
        width = 1
      ),
      strand = strand(gr),
      read_id = mcols(reads)$qname
    )
  
  return(five_prime_ends)
  
}

drop_unnecessary_seq <- function(x) {
  x <- keepSeqlevels(x, c("KF297339.1"))
  seqlengths(x) <- 237683
  return(x)
}

get_tsr <- function(x) {
  
  # Get coverage across genome of 5' ends 
  # The strand_coverage function is a wrapper around the GenomicRanges::coverage
  # function that allows for strand specific coverage 
  x.5primecov <- tsrDetectR::strand_coverage(x)
  x.5primecov.plus <- x.5primecov[[1]][[1]]
  x.5primecov.neg <- x.5primecov[[2]][[1]]
  
  # Call TSRs usings the `TSRFinder` based algorithm 
  maxtss.windows.plus <- tsrDetectR::tsr_maxtss(x.5primecov.plus, w = 51, threshold = 50)
  maxtss.windows.neg <- tsrDetectR::tsr_maxtss(x.5primecov.neg, w = 51, threshold = 50)
  
  # Create GRanges object to store TSRs 
  gr.plus <- GenomicRanges::GRanges(
    "KF297339.1", 
    maxtss.windows.plus, 
    strand="+",
    seqinfo = x@seqinfo)
  gr.neg <- GRanges(
    "KF297339.1", 
    maxtss.windows.neg, 
    strand="-",
    seqinfo = x@seqinfo)
  gr <- c(gr.plus, gr.neg)
  
  # Only keep the center position of each TSR, as that's the only part that's necessary 
  # for downstream analysis
  TSRs <- GenomicRanges::resize(gr, width = 1, fix = "center")
  
  return(TSRs)
}

sample_names <- c("mock", "4_hpi", "12_hpi", "24_hpi", "48_hpi", "72_hpi")
bam_files <- map(sample_names, ~list.files("results/alignments/", pattern = paste0("^", .x, ".+bam$"), full.names = T))

library(tsrDetectR)
five_prime_ends <- map(bam_files, extract_5prime_ends)
five_prime_ends <- map(five_prime_ends, drop_unnecessary_seq)

cmv_coverage <- map(five_prime_ends, ~coverage(.x)) |> 
  map(~as.vector(.x[[1]])) |> 
  map(sum)

tsrs <- map(five_prime_ends, get_tsr)
tsrs <- GRangesList(tsrs)
names(tsrs) <- sample_names

# Combine peaks into a single GRanges
tsrs_combined <- unlist(tsrs)
tsrs_combined$sample <- factor(tsrs_combined@ranges@NAMES, levels = sample_names)
tsrs_combined$id <- seq(length(tsrs_combined))
mcols(tsrs_combined) <- mcols(tsrs_combined) |> 
  data.frame() |> 
  pivot_wider(id_cols = id, names_from = sample, values_from = score) |> 
  select(all_of(sample_names))

tsrs_combined2 <- GenomicRanges::reduce(tsrs_combined, min.gapwidth = 10, with.revmap = T)

revmap <- mcols(tsrs_combined2)$revmap
mcols(tsrs_combined2) <- relist(tsrs_combined[unlist(revmap)], revmap) |> 
  map(~mcols(.x) |> data.frame() |> summarise_all(sum, na.rm=T)) |> 
  bind_rows() |> 
  select(mock, 
         `4hpi` = "X4_hpi", 
         `12hpi` = "X12_hpi",
         `24hpi` = "X24_hpi", 
         `48hpi` = "X48_hpi",
         `72hpi` = "X72_hpi")

# Remove RNA4.9 
tsrs_combined2 <- tsrs_combined2[!(
  tsrs_combined2@ranges@start == 94733 & 
    tsrs_combined2@strand == "+")]

tsrs_combined2 <- tsrs_combined2[!(
  tsrs_combined2@ranges@start == 94766 & 
    tsrs_combined2@strand == "-")]

tsr_scores <- mcols(tsrs_combined2) |> data.frame() |> tibble()
names(tsr_scores) <- sample_names
tsr_scores <- map2_df(tsr_scores, cmv_coverage, ~log(0.01 + 1e4*.x/.y))

library(ggalign)
ggheatmap(tsr_scores) +
  scale_fill_gradientn(colours = terrain.colors(5)) +
  anno_left() +
  align_dendro(k = 5, method = "ward.D")
