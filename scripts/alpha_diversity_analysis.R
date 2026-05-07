library(phyloseq)
library(tidyverse)
library(ggplot2)

# This script calculates alpha diversity metrics from a taxonomic
# count table and sample metadata.
#
# Expected input files:
#   data/species_counts.csv
#   data/metadata.txt
#
# Expected metadata columns:
#   Sample_ID
#   Group
#
# Expected count table:
#   First column = taxonomy / feature ID
#   Remaining columns = sample IDs matching metadata$Sample_ID
# ============================================================

# -----------------------------
# User-defined parameters
# -----------------------------

count_file <- "data/species_counts.csv"
metadata_file <- "data/metadata.txt"

sample_id_col <- "Sample_ID"
group_col <- "Group"

output_table <- "tables/alpha_diversity_summary.csv"
output_plot <- "figures/shannon_diversity.png"

# Optional: set to specific groups, e.g. c("Bran", "No food")
# Leave as NULL to include all groups
groups_to_keep <- NULL

# ============================================================
# Alpha Diversity Analysis
# ============================================================

# -----------------------------
# Load data
# -----------------------------

counts_raw <- read.csv(count_file, header = TRUE, check.names = FALSE)
metadata <- read.delim(metadata_file, header = TRUE, check.names = FALSE)

# Optional group filtering
if (!is.null(groups_to_keep)) {
  metadata <- metadata %>%
    filter(.data[[group_col]] %in% groups_to_keep)
}

samples_keep <- metadata[[sample_id_col]]

counts_filtered <- counts_raw %>%
  select(Feature = 1, all_of(samples_keep))

# -----------------------------
# Build taxonomy table
# -----------------------------

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

# -----------------------------
# Build phyloseq object
# -----------------------------

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

# -----------------------------
# Alpha diversity calculation
# -----------------------------

alpha_diversity <- estimate_richness(
  physeq,
  measures = c("Shannon", "Observed", "Simpson", "Chao")
)

alpha_results <- alpha_diversity %>%
  as.data.frame() %>%
  rownames_to_column(var = sample_id_col) %>%
  left_join(metadata, by = sample_id_col)

# -----------------------------
# Summary table
# -----------------------------

summary_table <- alpha_results %>%
  group_by(.data[[group_col]]) %>%
  summarise(
    mean_shannon = mean(Shannon, na.rm = TRUE),
    sd_shannon = sd(Shannon, na.rm = TRUE),
    mean_observed = mean(Observed, na.rm = TRUE),
    sd_observed = sd(Observed, na.rm = TRUE),
    mean_simpson = mean(Simpson, na.rm = TRUE),
    sd_simpson = sd(Simpson, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

write.csv(summary_table, output_table, row.names = FALSE)

# -----------------------------
# Visualisation
# -----------------------------

shannon_plot <- ggplot(
  alpha_results,
  aes(x = .data[[group_col]], y = Shannon, color = .data[[group_col]])
) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.8) +
  theme_bw() +
  labs(
    title = "Alpha Diversity by Group",
    x = NULL,
    y = "Shannon Index",
    color = "Group"
  )

ggsave(
  output_plot,
  plot = shannon_plot,
  width = 6,
  height = 5,
  dpi = 300
)

message("Alpha diversity analysis completed successfully.")
