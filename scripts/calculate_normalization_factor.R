library(dplyr)
library(tidyr)

# Read file
args <- commandArgs(TRUE)
read_depth_file = args[1]
mapped_reads <- read.delim(read_depth_file, sep = "\t")

# Reformat table, 1 row per sample
mapped_reads <- mapped_reads |> 
  separate(col = "Sample", 
           sep = "_(?=[A-Za-z0-9]+_extract$)|_(?=extract$)", # select last and second to last underscores, assuming each line ends with "extract"
           into = c("Sample", "Organism", "Extra")) |>
  select(-Extra) |> 
  pivot_wider(id_cols = "Sample", 
              names_from = "Organism", 
              values_from = "num_mapped_pairs")

# Rename last column 'spikein'
names(mapped_reads)[ncol(mapped_reads)] <- "spikein"

# Calculate the sum of mapped reads for each sample
mapped_reads <- mapped_reads |> 
  mutate(total = rowSums(across(-Sample)))

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
write.table(mapped_reads, file = read_depth_file, sep = "\t", row.names = F)