library(rtracklayer)
library(tsrDetectR)

files <- list.files("results/bed/", full.names = T)
peaks <- purrr::map(files, get_peaks)
peaks0 <- peaks

purrr::map2(
  peaks,
  list("mock", "4hpi", "12hpi", "24hpi", "48hpi", "72hpi"),
  ~export.bed(.x, paste0("results/peaks/", .y, ".bed")))
