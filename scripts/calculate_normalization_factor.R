library(dplyr)
library(tidyr)

# Read file
args <- commandArgs(TRUE)
mapped_reads <- read.delim(args[1], sep = "\t")

# Reformat table, 1 row per sample
mapped_reads <- mapped_reads |> 
  separate(col = "Sample", 
           sep = "_(?=[A-Za-z0-9]+_extract$)|_(?=extract$)", # select last and second to last underscores, assuming each line ends with "extract"
           into = c("Sample", "Organism", "Extra")) |>
  select(-Extra) |> 
  pivot_wider(id_cols = "Sample", 
              names_from = "Organism", 
              values_from = "num_mapped_pairs") |> 
  mutate(total = cmv + human + spikein)

# Calculate library Size Normalization
average_mapped_library_size = mean(mapped_reads$total)

mapped_reads <- mapped_reads |> 
  mutate(libsize_correction = average_mapped_library_size/total) |> 
  mutate(corrected_spikein = spikein*libsize_correction)

# Calculate spikein correction
average_corrected_spikein_size = mean(mapped_reads$corrected_spikein)

mapped_reads <- mapped_reads |> 
  mutate(spikein_correction = average_corrected_spikein_size/corrected_spikein) |> 
  mutate(total_correction = libsize_correction*spikein_correction)

# Save file
write.table(mapped_reads, file = args[1], sep = "\t", row.names = F)