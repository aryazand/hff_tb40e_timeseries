library(rtracklayer)

get_peaks <- function(bed_file) {

  x <- import.bed(bed_file) |> resize(width = 1, fix = "start")
  x.5primecov <- strand_coverage(x)
  x.5primecov.plus <- x.5primecov[[1]][[1]]
  x.5primecov.neg <- x.5primecov[[2]][[1]]


  maxtss.windows.plus <- tsr_maxtss(x.5primecov.plus, w = 31, threshold = 50)
  maxtss.windows.neg <- tsr_maxtss(x.5primecov.neg, w = 31, threshold = 50)

  gr.plus <- GRanges("KF297339.1", maxtss.windows.plus, strand="+",
                     seqinfo = x@seqinfo)
  gr.neg <- GRanges("KF297339.1", maxtss.windows.neg, strand="-",
                    seqinfo = x@seqinfo)

  x.peaks <- c(gr.plus, gr.neg)
  x.peaks <- x.peaks |> resize(width = 1, fix = "center")

  return(x.peaks)
}

file_pattern <- seq(52,63, 2) |> paste(collapse="|")
files <- list.files("~/hff_tb40_timeseries/results/bed/", full.names = T, pattern = file_pattern)

peaks <- purrr::map(files, get_peaks) |> GRangesList()
peaks0 <- peaks



purrr::map2(
  peaks,
  list("mock", "4hpi", "12hpi", "24hpi", "48hpi", "72hpi"),
  ~export.bed(.x, paste0("~/hff_tb40_timeseries/results/peaks/", .y, ".bed")))

purrr::map(c("mock", "4hpi", "12hpi", "24hpi", "48hpi", "72hpi"),
            ~paste0("~/hff_tb40_timeseries/results/peaks/", .x, ".bed")

export.bed(peaks[[1]], "~/test.bed")

#############
# Rename columns
#############

samples <- c("mock", "4hpi", "12hpi", "24hpi", "48hpi", "72hpi")
for(i in seq_along(peaks)) {
  names(mcols(peaks[[i]]))[1] <- samples[i]
}

peaks_combined <- do.call("c", peaks)
peaks_combined2 <- reduce(peaks_combined, min.gapwidth = 10, with.revmap = T)

revmap <- peaks_combined2$revmap
relist(peaks_combined[unlist(revmap)], revmap)


merge_peaks <- function(peaks1, peaks2, name1, name2, mingap = 16) {

  peaks.overlap <- c(peaks1, peaks2) |>
    reduce(min.gapwidth = mingap, with.revmap=TRUE)

  revmap <- mcols(peaks.overlap)$revmap

  relist(mcols(ir0)[unlist(revmap), ], revmap)

  return(peaks.overlap)

}

x <- GRanges()
for(i in 1:6) {
  x <- union(x, peaks[[i]])
}
