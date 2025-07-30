library(tidyverse)
library(VennDiagram)
library(RColorBrewer)

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# Load the pre-filtered significant gene set
sig_df <- read.csv("DA_results/sig_genes_combined.csv")

sig_df <- sig_df |> 
  mutate(
    adj_p <- p.adjust(pval, method = "BH")
  )

write.csv(sig_df, "DA_results/sig_genes_final.csv")

# Extract and standardize comparison labels
sig_df <- sig_df %>%
  mutate(
    comparison_raw = str_extract(combinations, "(?<=complete_uniref_).*"),
    comparison_label = str_split(comparison_raw, "vs", simplify = TRUE),
    comparison_label = paste(
      pmin(str_trim(comparison_label[, 1]), str_trim(comparison_label[, 2])),
      "vs",
      pmax(str_trim(comparison_label[, 1]), str_trim(comparison_label[, 2]))
    )
  )


# Define Venn groups
venn_panels <- list(
  B = c("Caries-CAa vs Healthy-CF", "Caries-CAi vs Healthy-CF", "Caries-CAs vs Healthy-CF"),
  C = c("Caries-CAa vs Healthy-CI", "Caries-CAi vs Healthy-CI", "Caries-CAs vs Healthy-CI"),
  D = c("Healthy-CF vs Healthy-CI", "Caries-CAa vs Healthy-CF", "Caries-CAa vs Healthy-CI"),
  E = c("Healthy-CF vs Healthy-CI", "Caries-CAi vs Healthy-CF", "Caries-CAi vs Healthy-CI"),
  F = c("Healthy-CF vs Healthy-CI", "Caries-CAs vs Healthy-CF", "Caries-CAs vs Healthy-CI"),
  G = c("Caries-CAa vs Caries-CAi", "Caries-CAa vs Caries-CAs", "Caries-CAi vs Caries-CAs"),
  H = c("Healthy-Baseline vs Healthy-Fl", "Healthy-Baseline vs Healthy-Fl-Ar", "Healthy-Fl vs Healthy-Fl-Ar"),
  I = c("Caries-Baseline vs Caries-Fl", "Caries-Baseline vs Caries-Fl-Ar", "Caries-Fl vs Caries-Fl-Ar"),
  J = c("Caries-Baseline vs Healthy-Baseline", "Caries-Fl vs Healthy-Fl", "Caries-Fl-Ar vs Healthy-Fl-Ar")
)

output_dir <- "DA_results/venn_panels/"
ensure_dir(output_dir)

all_sets_for_export <- list()
all_3way_intersections <- list()
# Generate Venn diagrams
for (panel in names(venn_panels)) {
  comparisons <- venn_panels[[panel]]
  
  sets <- lapply(comparisons, function(comp) {
    sig_df %>%
      filter(comparison_label == comp) %>%
      pull(Taxon) %>%
      unique()
  })
  names(sets) <- comparisons
  
  # Store for export
  all_sets_for_export <- c(all_sets_for_export, sets)
  common_all <- Reduce(intersect, sets)
  all_3way_intersections[[panel]] <- common_all
  if (all(sapply(sets, length) == 0)) {
    message("Skipping panel ", panel, " — all sets are empty.")
    next
  }
  
  insert_linebreaks <- function(label, words_per_line = 3) {
    words <- unlist(strsplit(label, " "))
    lines <- split(words, ceiling(seq_along(words) / words_per_line))
    paste(sapply(lines, paste, collapse = " "), collapse = "\n")
  }
  
  category_names_with_n <- mapply(
    function(name, taxa) {
      wrapped <- insert_linebreaks(name, words_per_line = 2)
      label_text <- paste0(wrapped, "\n(n = ", length(taxa), ")")
      label_text
    },
    names(sets), sets,
    SIMPLIFY = TRUE
  )
  
  file_path <- file.path(output_dir, paste0("Venn_Panel_", panel, ".png"))
  common_all <- Reduce(intersect, sets)
  message("Panel ", panel, ": 3-way intersection size = ", length(common_all))
  
  # Label angle and position
  if (panel == "J") {
    label_angle <- c(-10, -1, 10)
    label_dist <- c(0.15, 0.15, 0.05)
  } else {
    label_angle <- c(-22, 22, -180)
    label_dist <- c(0.2, 0.2, 0.2)
  }
  
  venn.diagram(
    x = sets,
    category.names = category_names_with_n,
    filename = file_path,
    imagetype = "png",
    height = 4000,
    width = 4000,
    resolution = 600,
    cex = 1.4,
    cat.cex = 1.4,
    fill = brewer.pal(3, "Set2"),
    cat.pos = label_angle,
    cat.dist = label_dist,
    margin = 0.1,
    col = "black",
    fontface = "bold",
    fontfamily = "Arial",
    cat.fontface = "bold",
    cat.fontfamily = "Arial",
    cat.col = rep("black", length(sets))
  )
  
  message("Saved: ", file_path)
}


max_len <- max(lengths(all_sets_for_export))
df_combined <- as.data.frame(lapply(all_sets_for_export, function(x) `length<-`(x, max_len)))
write.csv(df_combined, file.path(output_dir, "All_Venn_Sets_Combined.csv"), row.names = FALSE)

max_intersection_len <- max(lengths(all_3way_intersections))
df_intersections <- as.data.frame(lapply(all_3way_intersections, function(x) `length<-`(x, max_intersection_len)))
write.csv(df_intersections, file.path(output_dir, "All_Venn_Intersections_3way.csv"), row.names = FALSE)
