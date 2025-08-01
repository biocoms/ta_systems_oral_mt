# Load libraries
library(tidyverse)

# Set input paths
results_file <- "DA_results/sig_genes_combined.csv"
annotation_file <- "processed_data/result_annotation.csv" 

# Load data
results <- read_csv(results_file)
annotations <- read_csv(annotation_file)

# Clean condition labels if needed
results <- results %>%
  mutate(direction = case_when(lfc > 0 ~ "Up",
                               lfc < 0 ~ "Down",
                               TRUE ~ "NoChange"))


results <- results %>%
  mutate(
    comparison_raw = str_extract(combinations, "(?<=complete_uniref_).*"),
    comparison_label = str_split(comparison_raw, "vs", simplify = TRUE),
    condition_comparisons = paste(
      pmin(str_trim(comparison_label[, 1]), str_trim(comparison_label[, 2])),
      "vs",
      pmax(str_trim(comparison_label[, 1]), str_trim(comparison_label[, 2]))
    )
  )
# Merge with annotation
sig_annotated <- results %>%
  left_join(annotations, by = join_by("Taxon" == "gene_cluster"))

# Output 1: Summary table of number of DE genes per contrast
de_summary <- sig_annotated %>%
  group_by(condition_comparisons, direction) %>%
  summarise(n = n_distinct(Taxon), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  mutate(Total = Up + Down)

# View
print(de_summary)

# Output 2: List of annotated DE genes per contrast
de_genes_by_contrast <- sig_annotated %>%
  select(condition_comparisons, Taxon, Uniprot_Uniparc_ID, gene_id, protein_names_clean, direction) %>%
  arrange(condition_comparisons, direction, Taxon)

# Write output files
write_csv(de_summary, "DA_results/summary_DE_counts_by_contrast.csv")
write_csv(de_genes_by_contrast, "DA_results/annotated_DE_genes_by_contrast.csv")



# Define comparison categories
comparison_groups <- tribble(
  ~condition_comparisons,                            ~group,
  "Caries-CAa vs Healthy-CF",                       "Healthy_vs_Caries",
  "Caries-CAi vs Healthy-CF",                       "Healthy_vs_Caries",
  "Caries-CAs vs Healthy-CF",                       "Healthy_vs_Caries",
  "Caries-CAa vs Healthy-CI",                       "Healthy_vs_Caries",
  "Caries-CAi vs Healthy-CI",                       "Healthy_vs_Caries",
  "Caries-CAs vs Healthy-CI",                       "Healthy_vs_Caries",
  "Healthy-Baseline vs Healthy–Fl",                "Healthy_Intervention",
  "Healthy-Baseline vs Healthy-Fl-Ar",             "Healthy_Intervention",
  "Healthy-Fl vs Healthy-Fl-Ar",                   "Healthy_Intervention",
  "Caries-Baseline vs Caries–Fl",                  "Caries_Intervention",
  "Caries-Baseline vs Caries-Fl-Ar",               "Caries_Intervention",
  "Caries-Fl vs Caries-Fl-Ar",                     "Caries_Intervention",
  "Caries-Baseline vs Healthy-Baseline",           "Timepoint_Caries_vs_Healthy, Healthy_vs_Caries",
  "Caries-Fl vs Healthy-Fl",                       "Timepoint_Caries_vs_Healthy",
  "Caries-Fl-Ar vs Healthy-Fl-Ar", "Timepoint_Caries_vs_Healthy",
  "Caries-CAa vs Caries-CAi", "Spatial_Caries",
  "Caries-CAa vs Caries-CAs", "Spatial_Caries",
  "Caries-CAi vs Caries-CAs", "Spatial_Caries"
)



# Clean and extract meaningful TA family name
annotated_grouped <- sig_annotated %>%
  inner_join(comparison_groups, by = "condition_comparisons") %>%
  mutate(ta_family = protein_names_clean) %>%
  mutate(ta_family = if_else(ta_family == "", NA_character_, ta_family))

ta_family_patterns <- annotated_grouped %>%
  group_by(group, condition_comparisons, ta_family, direction) %>%
  summarise(n_genes = n_distinct(Taxon), .groups = "drop")

write_csv(ta_family_patterns, "DA_results/TA_family_patterns_by_group.csv")


# Optional: genes recurring across multiple comparisons
recurrence_table <- annotated_grouped %>%
  group_by(Taxon, protein_names_clean, ta_family) %>%
  summarise(n_contrasts = n_distinct(condition_comparisons),
            contrast_list = paste(unique(condition_comparisons), collapse = "; "),
            direction_list = paste(unique(direction), collapse = "; "),
            groups = paste(unique(group), collapse = "; "),
            .groups = "drop") %>%
  arrange(desc(n_contrasts))

write_csv(recurrence_table, "DA_results/Recurring_TA_clusters_across_contrasts.csv")


