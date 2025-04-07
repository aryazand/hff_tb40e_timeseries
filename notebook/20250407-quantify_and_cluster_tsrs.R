cmv_features <- rtracklayer::readGFF("data/KF297339.1.gff3")
cmv_features2 <- cmv_features |> data.frame() 
txdb <- txdbmaker::makeTxDbFromGFF("data/KF297339.1.gff3", format = "gff3")
cmv_promoters <- GenomicFeatures::promoters(txdb)

tsrs_combined4 <- tsrs_combined3
tsrs_combined4$nearest_promoter <- cmv_tx$tx_name[nearest(tsrs_combined4,cmv_promoters)]
tsrs_combined4$distance_to_promoter <- elementMetadata(distanceToNearest(tsrs_combined4, cmv_promoters))[["distance"]]

#########
# Identify Canonical Genes
##########

tsr_scores <- mcols(tsrs_combined4) |> data.frame() |> tibble()

IE.genes <- cmv_features2 |> 
  filter(grepl("IE[12]", product)) |> 
  pull(gene) |> 
  unique()

LTF.genes <- c("UL79", "UL87", "UL95")

tsr_scores <- tsr_scores %>% mutate(gene_category = 
  case_when(
    nearest_promoter %in% IE.genes & distance_to_promoter < 1000 ~ "Immediate Early",
    nearest_promoter %in% LTF.genes & distance_to_promoter < 1000 ~ "Late",
    .default = "other"
  )
)

names(tsr_scores) <- sub("^X", "", names(tsr_scores))
tsr_scores.normalized <- tsr_scores |> mutate_at(1:5, ~100*.x/sum(.x))
x <- tsr_scores.normalized |> mutate(tsr = as.character(row_number()))
x <- x |> pivot_longer(1:5, 
                       names_to = "time",
                       values_to = "score")
x$time <- factor(x$time, levels = sample_names[-1])

library(widyr)
y <- widely_kmeans(x, tsr, time, score, k = 5)
z <- left_join(x, y, by = "tsr")
z.normalized <- z |> group_by()



ggplot(z) + 
  geom_line(aes(time,score+0.000000001, group=tsr)) +
  facet_grid(gene_category~cluster) + 
  scale_y_log10()
