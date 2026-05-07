library(tidyverse)
library(KEGGREST)

# ============================================================
# KEGG Functional Profiling Analysis
# ============================================================

ko_file <- "data/ko_counts.txt"
metadata_file <- "data/metadata.txt"

sample_id_col <- "Sample_ID"
group_col <- "Group"

output_ko_counts <- "tables/kegg_id_counts_by_group.csv"
output_pathways <- "tables/kegg_pathway_summary.csv"

# -----------------------------
# Load data
# -----------------------------

ko_counts <- read.delim(ko_file, header = TRUE, check.names = FALSE)
metadata <- read.delim(metadata_file, header = TRUE, check.names = FALSE)

samples_keep <- metadata[[sample_id_col]]

ko_filtered <- ko_counts %>%
  select(KO = 1, all_of(samples_keep))

# -----------------------------
# KO presence by group
# -----------------------------

ko_long <- ko_filtered %>%
  pivot_longer(
    cols = -KO,
    names_to = sample_id_col,
    values_to = "Count"
  ) %>%
  left_join(metadata, by = sample_id_col)

total_unique_kos <- n_distinct(ko_long$KO)

ko_summary <- ko_long %>%
  group_by(KO, .data[[group_col]]) %>%
  summarise(
    total_count = sum(Count, na.rm = TRUE),
    present = total_count > 0,
    .groups = "drop"
  )

ko_presence_counts <- ko_summary %>%
  group_by(.data[[group_col]]) %>%
  summarise(
    KOs_present = sum(present),
    percentage_of_total = round((KOs_present / total_unique_kos) * 100, 2),
    .groups = "drop"
  )

write.csv(ko_presence_counts, output_ko_counts, row.names = FALSE)

# -----------------------------
# KEGG pathway mapping
# -----------------------------

ko_ids <- unique(gsub("^ko:", "", ko_summary$KO))

ko2path_raw <- keggLink("pathway", "ko")

ko2path <- tibble(
  KO = gsub("^ko:", "", names(ko2path_raw)),
  Pathway = gsub("^path:", "", as.character(ko2path_raw))
) %>%
  filter(KO %in% ko_ids)

pathway_info_raw <- keggList("pathway")

pathway_info <- tibble(
  Pathway = gsub("^path:", "", names(pathway_info_raw)),
  Description = as.character(pathway_info_raw)
)

ko_pathway <- ko_summary %>%
  mutate(KO_clean = gsub("^ko:", "", KO)) %>%
  left_join(ko2path, by = c("KO_clean" = "KO")) %>%
  left_join(pathway_info, by = "Pathway")

# -----------------------------
# Core and unique functional pathways
# -----------------------------

pathway_presence <- ko_pathway %>%
  filter(!is.na(Description)) %>%
  group_by(Description, .data[[group_col]]) %>%
  summarise(
    present = any(present),
    n_KOs = n_distinct(KO[present]),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = .data[[group_col]],
    values_from = present,
    values_fill = FALSE
  )

group_columns <- setdiff(colnames(pathway_presence), "Description")

pathway_summary <- pathway_presence %>%
  rowwise() %>%
  mutate(
    groups_present = paste(group_columns[c_across(all_of(group_columns))], collapse = ";"),
    n_groups_present = sum(c_across(all_of(group_columns))),
    Category = case_when(
      n_groups_present == length(group_columns) ~ "Core",
      n_groups_present == 1 ~ paste0("Unique_to_", groups_present),
      TRUE ~ "Partially_shared"
    )
  ) %>%
  ungroup()

write.csv(pathway_summary, output_pathways, row.names = FALSE)

message("KEGG functional profiling analysis completed successfully.")
