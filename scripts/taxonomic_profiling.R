library(phyloseq)
library(tidyverse)
library(ggplot2)

# ============================================================
#  Taxonomic Profiling Analysis
# ============================================================

count_file <- "data/species_counts.csv"
metadata_file <- "data/metadata.txt"

sample_id_col <- "Sample_ID"
group_col <- "Group"

output_table <- "tables/dominant_taxa_by_group.csv"
output_plot <- "figures/top_taxa_relative_abundance.png"

top_n <- 10

counts_raw <- read.csv(count_file, header = TRUE, check.names = FALSE)
metadata <- read.delim(metadata_file, header = TRUE, check.names = FALSE)

samples_keep <- metadata[[sample_id_col]]

counts_filtered <- counts_raw %>%
  select(Feature = 1, all_of(samples_keep))

otu_df <- data.frame(counts_filtered)
colnames(otu_df)[1] <- "Feature"

taxonomy_df <- data.frame(Feature = otu_df$Feature)

taxonomy_split <- taxonomy_df %>%
  separate(
    col = Feature,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = ";",
    remove = FALSE,
    fill = "right",
    extra = "merge"
  )

otu_matrix <- otu_df %>%
  column_to_rownames("Feature") %>%
  as.matrix()

taxonomy_matrix <- taxonomy_split %>%
  column_to_rownames("Feature") %>%
  as.matrix()

metadata_df <- metadata %>%
  column_to_rownames(sample_id_col)

physeq <- phyloseq(
  otu_table(otu_matrix, taxa_are_rows = TRUE),
  tax_table(taxonomy_matrix),
  sample_data(metadata_df)
)

# Convert to relative abundance
physeq_relative <- transform_sample_counts(
  physeq,
  function(x) x / sum(x)
)

relative_df <- psmelt(physeq_relative)

# Create taxonomic string
relative_df <- relative_df %>%
  mutate(
    Taxonomic_string = paste(
      Kingdom, Phylum, Class, Order, Family, Genus, Species,
      sep = ";"
    )
  )

# Dominant taxon per group
dominant_taxa <- relative_df %>%
  group_by(.data[[group_col]], Taxonomic_string) %>%
  summarise(
    mean_relative_abundance = mean(Abundance) * 100,
    .groups = "drop"
  ) %>%
  group_by(.data[[group_col]]) %>%
  slice_max(mean_relative_abundance, n = 1) %>%
  ungroup()

write.csv(dominant_taxa, output_table, row.names = FALSE)

# Top taxa overall
top_taxa <- relative_df %>%
  group_by(Taxonomic_string) %>%
  summarise(
    mean_relative_abundance = mean(Abundance),
    .groups = "drop"
  ) %>%
  slice_max(mean_relative_abundance, n = top_n) %>%
  pull(Taxonomic_string)

plot_df <- relative_df %>%
  mutate(
    Taxon_grouped = if_else(
      Taxonomic_string %in% top_taxa,
      Taxonomic_string,
      "Other"
    )
  ) %>%
  group_by(.data[[group_col]], Taxon_grouped) %>%
  summarise(
    mean_relative_abundance = mean(Abundance) * 100,
    .groups = "drop"
  )

taxa_plot <- ggplot(
  plot_df,
  aes(
    x = .data[[group_col]],
    y = mean_relative_abundance,
    fill = Taxon_grouped
  )
) +
  geom_col() +
  theme_bw() +
  labs(
    title = "Top Taxa Relative Abundance",
    x = NULL,
    y = "Mean Relative Abundance (%)",
    fill = "Taxon"
  )

ggsave(
  output_plot,
  plot = taxa_plot,
  width = 9,
  height = 6,
  dpi = 300
)

message("Taxonomic profiling analysis completed successfully.")
