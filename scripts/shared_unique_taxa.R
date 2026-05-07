library(phyloseq)
library(tidyverse)

# ============================================================
# Shared and Unique Taxa Analysis
# ============================================================

count_file <- "data/species_counts.csv"
metadata_file <- "data/metadata.txt"

sample_id_col <- "Sample_ID"
group_col <- "Group"

output_file <- "tables/shared_unique_taxa.csv"

counts_raw <- read.csv(count_file, header = TRUE, check.names = FALSE)
metadata <- read.delim(metadata_file, header = TRUE, check.names = FALSE)

samples_keep <- metadata[[sample_id_col]]

counts_filtered <- counts_raw %>%
  select(Feature = 1, all_of(samples_keep))

long_counts <- counts_filtered %>%
  pivot_longer(
    cols = -Feature,
    names_to = sample_id_col,
    values_to = "Count"
  ) %>%
  left_join(metadata, by = sample_id_col)

taxa_presence <- long_counts %>%
  group_by(Feature, .data[[group_col]]) %>%
  summarise(
    total_count = sum(Count, na.rm = TRUE),
    present = total_count > 0,
    .groups = "drop"
  ) %>%
  filter(present)

presence_wide <- taxa_presence %>%
  mutate(value = TRUE) %>%
  select(Feature, .data[[group_col]], value) %>%
  pivot_wider(
    names_from = .data[[group_col]],
    values_from = value,
    values_fill = FALSE
  )

group_columns <- setdiff(colnames(presence_wide), "Feature")

shared_unique <- presence_wide %>%
  rowwise() %>%
  mutate(
    groups_present = paste(group_columns[c_across(all_of(group_columns))], collapse = ";"),
    n_groups_present = sum(c_across(all_of(group_columns)))
  ) %>%
  ungroup() %>%
  mutate(
    Category = case_when(
      n_groups_present == length(group_columns) ~ "Shared",
      n_groups_present == 1 ~ paste0("Unique_to_", groups_present),
      TRUE ~ "Partially_shared"
    )
  )

write.csv(shared_unique, output_file, row.names = FALSE)

message("Shared and unique taxa analysis completed successfully.")
