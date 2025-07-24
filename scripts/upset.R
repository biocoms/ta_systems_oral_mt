# Load required libraries
library(UpSetR)
library(dplyr)
library(tidyr)
library(grid)

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# Define valid conditions with spaces around "vs"
valid_dieguez_conditions <- c(
  "Healthy-Baseline vs Caries-Baseline", "Healthy-Fl vs Caries-Fl", "Healthy-Fl-Ar vs Caries-Fl-Ar",
  "Healthy-Baseline vs Healthy-Fl", "Healthy-Fl vs Healthy-Fl-Ar", "Healthy-Baseline vs Healthy-Fl-Ar",
  "Caries-Baseline vs Caries-Fl", "Caries-Fl vs Caries-Fl-Ar", "Caries-Baseline vs Caries-Fl-Ar"
)

valid_ev_conditions <- c(
  "Healthy-CF vs Healthy-CI", "Healthy-CF vs Caries-CAa", "Healthy-CF vs Caries-CAi", "Healthy-CF vs Caries-CAs",
  "Healthy-CI vs Caries-CAa", "Healthy-CI vs Caries-CAi", "Healthy-CI vs Caries-CAs",
  "Caries-CAa vs Caries-CAi", "Caries-CAi vs Caries-CAs", "Caries-CAa vs Caries-CAs"
)

valid_conditions <- c(valid_dieguez_conditions, valid_ev_conditions)

tadb3_ta <- read_tsv("tadb3/tadb3_tox_antitox_main_table.tsv")

tadb3_toxins_sig <- tadb3_ta %>%
  filter(Toxin_Antitoxin == "Toxin") %>%
  pull(UniRef90_ID) %>%
  unique()

tadb3_antitoxins_sig <- tadb3_ta %>%
  filter(Toxin_Antitoxin == "Antitoxin") %>%
  pull(UniRef90_ID) %>%
  unique()

sig_genes_binary <- sig_genes_combined %>%
  mutate(condition_only = sub("^[^_]+_[^_]+_[^_]+_", "", combinations)) %>%
  mutate(condition_only = gsub("vs", " vs ", condition_only)) %>%
  mutate(
    standard_condition = ifelse(
      condition_only %in% valid_conditions, 
      condition_only,
      ifelse(
        gsub("^(.*) vs (.*)$", "\\2 vs \\1", condition_only) %in% valid_conditions,
        gsub("^(.*) vs (.*)$", "\\2 vs \\1", condition_only),
        NA_character_
      )
    )
  ) %>%
  filter(!is.na(standard_condition))

print(sig_genes_binary)

presence_matrix <- sig_genes_binary %>%
  select(Taxon, standard_condition) %>%
  distinct() %>%
  pivot_wider(names_from = standard_condition, values_from = standard_condition, values_fill = list(standard_condition = 0),
              values_fn = length) %>%
  mutate(across(-Taxon, ~ifelse(. > 0, 1, 0)))

presence_matrix_binary <- presence_matrix %>%
  mutate(across(-Taxon, ~ ifelse(. > 0, 1, 0)))  # Convert to 1/0 binary

binary_matrix <- presence_matrix_binary[, -c(1, ncol(presence_matrix_binary))]  # Exclude Taxon & Gene
binary_matrix <- apply(binary_matrix, 2, as.integer)  # Ensure numeric
binary_matrix <- as.data.frame(binary_matrix)  # Convert back to dataframe

# Ensure column names are properly set
colnames(binary_matrix) <- colnames(presence_matrix_binary)[-c(1, ncol(presence_matrix_binary))]
binary_matrix <- cbind(Taxon = presence_matrix_binary$Taxon, binary_matrix)

set.metadata <- data.frame(
  sets = colnames(binary_matrix)[-1],  # Extract comparison names
  category = ifelse(colnames(binary_matrix)[-1] %in% valid_dieguez_conditions, "dieguez",
                    ifelse(colnames(binary_matrix)[-1] %in% valid_ev_conditions, "ev", "Other"))
)
color_map <- c("dieguez" = "#E82561", "ev" = "purple", "Other" = "gray")
ensure_dir("DA_results/upset_plot")


png("DA_results/upset_plot/upset_plot.png", width = 5000, height = 3500, res = 600)

print(
  UpSetR::upset(
    binary_matrix,
    sets = colnames(binary_matrix)[-1],
    keep.order = TRUE,  
    order.by = "freq",  
    main.bar.color = "royalblue3",  
    sets.bar.color = color_map[set.metadata$category],  
    matrix.color = "#B82132",  
    text.scale = c(1.2, 1.5, 1.2, 0.8, 0.8, 1.2),  
    mainbar.y.label = "Gene Cluster Sets Intersections",  
    sets.x.label = "Significant genes \n in each comparison",
    set_size.show = TRUE,  
    set_size.scale_max = 100,
    queries = list(
      list(query = elements, params = list("Taxon", tadb3_toxins_sig), color = "#E52020", active = TRUE, query.name = "Toxins"),
      list(query = elements, params = list("Taxon", tadb3_antitoxins_sig), color = "#00FF9C", active = TRUE, query.name = "Antitoxins")
    ),
    query.legend = "bottom"
  )
)

dev.off()

#print(upset_plot)
