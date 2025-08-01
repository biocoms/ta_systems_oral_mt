# # Install packages
# install.packages(c(
#   "tidyr", "tidyverse", "dplyr", "VennDiagram", "pheatmap", 
#   "ggplot2", "patchwork", "gridExtra", "RColorBrewer", "UpSetR"
# ))
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("phyloseq", "ANCOMBC", "ComplexHeatmap", "circlize", "microbiome"))

# Load necessary libraries
library(tidyr)
library(tidyverse)
library(dplyr)
library(VennDiagram)
library(phyloseq)
library(ANCOMBC)
library(pheatmap)
library(ggplot2)
library(ComplexHeatmap)
library(patchwork)
library(gridExtra)  
library(circlize)
library(RColorBrewer)
library(UpSetR)
library(microbiome)

# Load datasets
dieguez_gene <- read.csv("raw_data/dieguez_genefamilies.tsv", sep="\t")
ev_gene <- read.csv("raw_data/ev_genefamilies.tsv", sep = "\t")
uniref <- read.csv("raw_data/uniref90_toxin_toxic.tsv", sep ="\t")
metadata_dieguez <- read_csv("raw_data/metadata_dieguez.csv")
metadata_ev <- read_csv("raw_data/metadata_ev.csv")
metadata_dieguez <- as.data.frame(metadata_dieguez)
metadata_ev <- as.data.frame(metadata_ev)

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# Data Cleaning Gene Data
# Preprocess dieguez_gene
dieguez_gene <- dieguez_gene %>%
  filter(rowSums(select_if(., is.numeric)) > 0)
colnames(dieguez_gene)[1] <- "gene_family"
colnames(dieguez_gene) <- gsub("_unaligned_Abundance.RPKs", "", colnames(dieguez_gene))
dieguez_gene <- dieguez_gene %>%
  filter(gene_family != "UNMAPPED")

# Preprocess ev_gene
colnames(ev_gene)[1] <- "gene_family"
colnames(ev_gene) <- gsub("_unaligned_Abundance.RPKs", "", colnames(ev_gene))
ev_gene <- ev_gene %>%
  filter(gene_family != "UNMAPPED")

dim(dieguez_gene)
dim(ev_gene)
head(dieguez_gene, 10)
head(ev_gene, 10)

# For dieguez_gene dataset
filtered_dieguez_gene <- subset(dieguez_gene, grepl("^UniRef90_[^|]+$", gene_family))
filtered_dieguez_gene_org <- subset(dieguez_gene, grepl("^UniRef90_[^|]+\\|", gene_family))

# For ev_gene dataset
filtered_ev_gene <- subset(ev_gene, grepl("^UniRef90_[^|]+$", gene_family))
filtered_ev_gene_org <- subset(ev_gene, grepl("^UniRef90_[^|]+\\|", gene_family))

# Displaying the results
dim(filtered_dieguez_gene)
head(filtered_dieguez_gene)
dim(filtered_ev_gene)
head(filtered_ev_gene)
dim(filtered_dieguez_gene_org)
head(filtered_dieguez_gene_org)
dim(filtered_ev_gene_org)
head(filtered_ev_gene_org)

# Remove duplicate gene families in dieguez and ev datasets
dieguez_gene_samples_unique <- filtered_dieguez_gene %>%
  distinct(gene_family, .keep_all = TRUE)  # Keep only the first occurrence
ev_gene_samples_unique <- filtered_ev_gene %>%
  distinct(gene_family, .keep_all = TRUE)  # Keep only the first occurrence
rownames(dieguez_gene_samples_unique)<- dieguez_gene_samples_unique$gene_family
dieguez_counts <- as.matrix(dieguez_gene_samples_unique[, -1]) 
rownames(dieguez_counts)<- dieguez_gene_samples_unique$gene_family
rownames(ev_gene_samples_unique)<- ev_gene_samples_unique$gene_family
ev_counts <- as.matrix(ev_gene_samples_unique[, -1])
rownames(ev_counts)<- ev_gene_samples_unique$gene_family


# Save the processed files
ensure_dir("processed_data")
write.csv(dieguez_counts, "processed_data/processed_dieguez_genes.csv", row.names = TRUE)
write.csv(ev_counts, "processed_data/processed_ev_genes.csv", row.names = TRUE)
write.csv(filtered_dieguez_gene_org, "processed_data/processed_dieguez_gene_org.csv", row.names = FALSE)
write.csv(filtered_ev_gene_org, "processed_data/processed_ev_gene_org.csv", row.names = FALSE)


# Overlaps with UniRef90 Toxin-related clusters
# Extract unique UniRef90 IDs for overlap comparison
dieguez_uniref_ids <- unique(rownames(dieguez_counts))
ev_uniref_ids <- unique(rownames(ev_counts))
uniref_ids <- unique(uniref$`Cluster.ID`)

# Find overlaps
overlap_dieguez_uniref <- intersect(dieguez_uniref_ids, uniref_ids)
cat("Number of overlapping Ids between dieguez and UniRef90 Toxins list:",length(overlap_dieguez_uniref), "\n")
overlap_ev_uniref <- intersect(ev_uniref_ids, uniref_ids)
cat("Number of overlapping Ids between ev and UniRef90 Toxins list:",length(overlap_ev_uniref), "\n")

# Find union and intersection of overlaps
union_uniref_ids <- union(overlap_dieguez_uniref, overlap_ev_uniref)
intersection_dieguez_ev <- intersect(overlap_dieguez_uniref, overlap_ev_uniref)
cat("Number of overlapping UniRef90 IDs between dieguez and ev:", length(intersection_dieguez_ev), "\n")

# Filter datasets based on overlaps
# For dieguez - overlap with UniRef but not overlapping with ev
dieguez_not_in_ev <- dieguez_gene_samples_unique %>%
  filter(gene_family %in% overlap_dieguez_uniref & !(gene_family %in% intersection_dieguez_ev))
cat("Number of UniRef90 Ids overlapping between dieguez and UniRef90 Toxins list but not in ev:",length(dieguez_not_in_ev$gene_family), "\n")

# For ev - overlap with UniRef but not overlapping with dieguez
ev_not_in_dieguez <- ev_gene_samples_unique %>%
  filter(gene_family %in% overlap_ev_uniref & !(gene_family %in% intersection_dieguez_ev))
cat("Number of UniRef90 Ids overlapping between ev and UniRef90 Toxins list but not in dieguez:",length(ev_not_in_dieguez$gene_family), "\n")

# For Intersection - overlap between both dieguez and ev with UniRef
intersection_dieguez <- dieguez_gene_samples_unique %>%
  filter(gene_family %in% intersection_dieguez_ev)

intersection_ev <- ev_gene_samples_unique %>%
  filter(gene_family %in% intersection_dieguez_ev)

# Save the filtered datasets
ensure_dir("uniref_tox_abundance")
write.csv(dieguez_not_in_ev, "uniref_tox_abundance/dieguez_final_uniref90_mapped.csv", row.names = FALSE)
write.csv(ev_not_in_dieguez, "uniref_tox_abundance/ev_final_uniref90_mapped.csv", row.names = FALSE)
write.csv(intersection_dieguez, "uniref_tox_abundance/intersection_dieguez_final_uniref90_mapped.csv", row.names = FALSE)
write.csv(intersection_ev, "uniref_tox_abundance/intersection_ev_final_uniref90_mapped.csv", row.names = FALSE)


# Venn Diagrams of the overlaps
# Prepare data for Venn Diagram
venn_data <- list(
  dieguez = overlap_dieguez_uniref,
  ev = overlap_ev_uniref
)

# Generate the Venn Diagram with improved colors
venn_plot <- venn.diagram(
  x = venn_data,
  category.names = c("Dieguez & UniRef90 TAs", "Ev & UniRef90 TAs"),
  filename = NULL,  # NULL to draw in R
  col = "black",    # Border color
  fill = c("skyblue", "pink1"),
  alpha = 0.5,      # Transparency level
  cex = 1,          # Font size for numbers
  fontface = "bold", # Bold text for numbers
  fontfamily = "sans", # Font family for numbers
  cat.cex = 0.7,      # Font size for category labels
  cat.fontface = "bold", # Bold text for categories
  cat.fontfamily = "sans", # Font family for categories
  cat.col = "black", # Matching darker colors for category labels
  cat.dist = c(0.1, 0.1), # Adjust distances of labels from circles
  margin = 0.15    # Margins around the diagram
)

# Save the Venn diagram to a file
ensure_dir("venn_diagrams")
png("venn_diagrams/venn_dieguez_ev_uniref.png", width = 2500, height = 2500, res = 600)
grid.draw(venn_plot)
dev.off()

# Prepare subset for dieguez
dieguez_complete_uniref <- rbind(dieguez_not_in_ev, intersection_dieguez)
dieguez_complete_uniref_counts <- as.matrix(dieguez_complete_uniref[, -1])  # Exclude gene_family column
rownames(dieguez_complete_uniref_counts) <- dieguez_complete_uniref$gene_family

# Prepare subset for ev
ev_complete_uniref <- rbind(ev_not_in_dieguez, intersection_ev)
ev_complete_uniref_counts <- as.matrix(ev_complete_uniref[, -1])  # Exclude gene_family column
rownames(ev_complete_uniref_counts) <- ev_complete_uniref$gene_family

# Ensure metadata alignment
rownames(metadata_ev) <- metadata_ev$run
rownames(metadata_dieguez) <- metadata_dieguez$run

# Phyloseq for dieguez
physeq_dieguez <- phyloseq(
  otu_table(dieguez_complete_uniref_counts, taxa_are_rows = TRUE),
  sample_data(metadata_dieguez)
)

# Phyloseq for ev
physeq_ev <- phyloseq(
  otu_table(ev_complete_uniref_counts, taxa_are_rows = TRUE),
  sample_data(metadata_ev)
)

standardize_filename <- function(name) {
  name |>
    gsub("_complete_uniref|_uniref", "", x = _) |>  # remove suffixes
    gsub("\\.csv$", "", x = _) |>                  # remove .csv extension
    gsub("vs", "_vs_", x = _) |>                   # standardize 'vs'
    gsub("-", "_", x = _) |>                       # replace dash with underscore
    gsub("\\s+", "_", x = _) |>                    # collapse spaces to single underscore
    gsub("_+", "_", x = _) |>                      # collapse multiple underscores
    gsub("^_|_$", "", x = _)                       # trim leading/trailing _
}


make_plot_title <- function(file_path) {
  title <- basename(file_path) |>
    gsub("_complete_uniref|_uniref", "", x = _) |>   # strip suffixes
    gsub("\\.csv$", "", x = _) |>                   # strip .csv
    gsub("_vs_", " vs ", x = _) |>                  # format 'vs'
    gsub("_", " ", x = _) |>                        # underscores to spaces
    gsub("(?<=[a-z])(?=[A-Z])", " ", x = _, perl=TRUE) |> # spacing CamelCase (if any)
    tools::toTitleCase()                            # proper title case
  return(title)
}

# DA analysis
run_ancombc_analysis <- function(datasets, metadata, output_dir = "DA_results/") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Create subdirectories for complete and significant results
  complete_dir <- file.path(output_dir, "da_complete")
  sig_dir <- file.path(output_dir, "da_sig")
  if (!dir.exists(complete_dir)) dir.create(complete_dir, recursive = TRUE)
  if (!dir.exists(sig_dir)) dir.create(sig_dir, recursive = TRUE)
  
  sig_genes_list <- list()
  
  for (dataset in datasets) {
    dataset_name <- dataset$name
    physeq <- dataset$physeq
    
    # Extract unique conditions
    conditions <- unique(metadata$condition)
    
    # Generate all pairwise combinations
    comb_pairs <- combn(conditions, 2, simplify = FALSE)
    
    for (comb in comb_pairs) {
      condition1 <- comb[1]
      condition2 <- comb[2]
      
      # Filter metadata for the two conditions
      meta_filtered <- metadata %>% filter(condition %in% c(condition1, condition2))
      
      # Check if there are samples for these conditions
      if (nrow(meta_filtered) == 0) {
        cat("No samples found for conditions:", condition1, "and", condition2, "in dataset:", dataset_name, "\n")
        next
      }
      
      # Update metadata for the subset
      meta_filtered$condition <- factor(meta_filtered$condition, levels = c(condition1, condition2))
      
      # Subset the phyloseq object
      physeq_subset <- prune_samples(sample_names(physeq) %in% meta_filtered$run, physeq)
      
      # Run ANCOMBC
      ancombc_res <- tryCatch({
        ancombc(
          data = physeq_subset,
          formula = "condition"
        )
      }, error = function(e) {
        cat("Error in ANCOMBC for conditions:", condition1, "and", condition2, "in dataset:", dataset_name, "-", e$message, "\n")
        NULL
      })
      
      # Skip if ANCOMBC fails
      if (is.null(ancombc_res)) next
      
      # Extract results
      res <- ancombc_res$res
      condition_col <- colnames(res$p_val)[3]  # Extract the correct condition column dynamically
      
      # Create result dataframe
      result_df <- data.frame(
        Taxon = res$p_val$taxon,  # Extract taxon names
        pval = res$p_val[[condition_col]],  # Extract p-values for the condition
        qval = res$q_val[[condition_col]],  # Extract q-values for the condition
        lfc = res$lfc[[condition_col]]      # Extract log-fold change for the condition
      )
      
      # Clean combination name
      combo_raw <- paste(condition1, "vs", condition2)
      combo_clean <- gsub("\\s+", "", combo_raw)
      filename_clean <- paste0(dataset_name, "_", combo_clean)
      
      # Save files
      full_results_file <- file.path(complete_dir, paste0(standardize_filename(filename_clean), ".csv"))
      
      write.csv(result_df, file = full_results_file, row.names = FALSE)
      
      # Filter significant taxa (q-value or p-value < 0.05)
      significant_df <- result_df %>% filter(qval < 0.05 | pval < 0.05)
      
      # Save significant results
      sig_results_file <- file.path(sig_dir, paste0(standardize_filename(filename_clean), ".csv"))
      write.csv(significant_df, file = sig_results_file, row.names = FALSE)
      
      # Store significant results in a list for final dataframe
      if (nrow(significant_df) > 0) {
        significant_df <- significant_df %>%
          mutate(combinations = paste0(dataset_name, "_", gsub(" ", "", condition1), "vs", gsub(" ", "", condition2))) %>%
          select(Taxon, combinations, lfc, pval, qval)
        
        sig_genes_list <- append(sig_genes_list, list(significant_df))
      }
    }
  }
  
  # Convert list to a single dataframe
  sig_genes_df <- do.call(rbind, sig_genes_list)
  
  # Return the combined significant genes dataframe
  return(sig_genes_df)
}

# Prepare datasets for dieguez and ev
datasets_dieguez <- list(
  list(name = "dieguez_complete_uniref", physeq = physeq_dieguez)
)

datasets_ev <- list(
  list(name = "ev_complete_uniref", physeq = physeq_ev)
)

# Run for both dieguez and ev and combine results
sig_genes_dieguez <- run_ancombc_analysis(datasets_dieguez, metadata_dieguez)
sig_genes_ev <- run_ancombc_analysis(datasets_ev, metadata_ev)

# Combine results from both analyses
sig_genes_combined <- rbind(sig_genes_dieguez, sig_genes_ev)

# Display final significant results
print(sig_genes_combined)

# Save combined significant genes dataframe
write.csv(sig_genes_combined, file = "DA_results/sig_genes_combined.csv", row.names = FALSE)

# Volcano plots for the complete datasets
create_volcano_plot <- function(file_path, output_name, pval_cutoff = 0.05, lfc_cutoff = 1) {
  # Read the data
  data <- read.csv(file_path)
  
  # Ensure required columns exist
  if (!all(c("Taxon", "pval", "qval", "lfc") %in% colnames(data))) {
    cat("Skipping:", file_path, "- Required columns not found.\n")
    return(NULL)
  }
  
  # Create log10-transformed p-values
  data$log10_pval <- -log10(data$pval)
  
  # Define significance based on p-value & LFC
  data$significant <- ifelse(data$pval < pval_cutoff, "Significant", "Not Significant")
  
  # Define highly significant genes (optional)
  #data$top_genes <- ifelse(data$significant == "Significant" & data$pval < (pval_cutoff / 10), as.character(data$Taxon), NA)
  wrap_by_words <- function(text, words_per_line = 5) {
    words <- unlist(strsplit(text, "\\s+"))
    lines <- split(words, ceiling(seq_along(words) / words_per_line))
    wrapped_text <- sapply(lines, paste, collapse = " ")
    paste(wrapped_text, collapse = "\n")
  }
  
  # Create the volcano plot
  p <- ggplot(data, aes(x = lfc, y = log10_pval, color = significant)) +
    geom_point(alpha = 1.0, size = 4) +
    scale_color_manual(values = c("gray", "red"), name = "Significance") +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dotted", color = "blue") +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = "dotted", color = "blue") +
    theme_minimal() +
    labs(
      title = paste(wrap_by_words(make_plot_title(file_path), 7)),
      x = "Log Fold Change",
      y = "-log10(p-value)"
    ) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 12)
    )
  
  # Save plots in TIFF and PNG formats
  # ggsave(paste0(output_name, ".tiff"), plot = p, dpi = 600, width = 5, height = 5, units = "in")
  ggsave(paste0(output_name, ".png"), plot = p, dpi = 600, width = 5, height = 5, units = "in")
  
  cat("Volcano plot saved for:", file_path, "\n")
}

# Process all relevant files
all_files <- list.files("DA_results/da_complete", pattern = ".csv$", full.names = TRUE)

# Loop through each file and create plots
for (file_path in all_files) {
  ensure_dir("DA_results/volcano_plot")
  # Generate an output name based on the file name
  output_name <- paste0("DA_results/volcano_plot/",
                        standardize_filename(basename(file_path)))

  # Create the volcano plot
  create_volcano_plot(file_path, output_name)
}


# Function to generate heatmaps from DA files
generate_all_heatmaps <- function(da_files, counts_matrix, metadata, dataset_name, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  metadata$condition <- metadata$condition %>%
    trimws() %>%
    stringr::str_replace_all("\\s{2,}", " ")
  
  metadata_conditions <- unique(metadata$condition)
  
  clean_condition <- function(cond_string) {
    cond_string %>%
      gsub("_", " - ", .) %>%
      gsub("\\s{2,}", " ", .) %>%
      stringr::str_trim() %>%
      stringr::str_replace_all("\\bFl Ar\\b", "Fl-Ar") %>%
      stringr::str_replace_all("\\bCI\\b", "CI") %>%
      stringr::str_replace_all("\\bCF\\b", "CF") %>%
      tools::toTitleCase()
  }
  
  
  for (da_file in da_files) {
    file_base <- basename(da_file)
    comparison <- file_base %>%
      gsub("^dieguez_|^ev_", "", .) %>%
      gsub("\\.csv$", "", .)
    
    cat("\nProcessing file:", file_base, "-> Comparison:", comparison, "\n")
    
    da_data <- readr::read_csv(da_file, show_col_types = FALSE)
    significant_taxa <- unique(da_data$Taxon)
    sig_counts <- counts_matrix[rownames(counts_matrix) %in% significant_taxa, , drop = FALSE]
    
    if (nrow(sig_counts) == 0) {
      cat("No significant genes found for", comparison, "- Skipping heatmap\n")
      next
    }
    
    sig_counts_log <- log2(sig_counts + 1)
    sig_counts_scaled <- t(scale(t(sig_counts_log)))
    
    conditions <- unlist(strsplit(comparison, "_vs_"))
    conditions <- lapply(conditions, clean_condition)
    conditions <- unlist(conditions)
    
    if (length(conditions) != 2) {
      cat("Invalid comparison format (expected 2 conditions):", comparison, "\n")
      next
    }
    
    cat("Extracted Conditions:", conditions[1], "vs", conditions[2], "\n")
    cat("Available metadata conditions:\n")
    print(metadata_conditions)
    
    match_condition <- function(cond) {
      exact <- metadata_conditions[metadata_conditions == cond]
      if (length(exact) > 0) return(exact[1])
      fuzzy <- metadata_conditions[agrep(cond, metadata_conditions, ignore.case = TRUE, max.distance = 0.2)]
      if (length(fuzzy) > 0) return(fuzzy[1])
      return(NA)
    }
    
    condition1 <- match_condition(conditions[1])
    condition2 <- match_condition(conditions[2])
    
    if (is.na(condition1) || is.na(condition2) || condition1 == condition2) {
      cat("Could not match both conditions or they are identical:", condition1, condition2, "\n")
      next
    }
    
    cat("Final Matched Conditions:", condition1, "vs", condition2, "\n")
    
    meta_filtered <- metadata %>%
      dplyr::filter(condition %in% c(condition1, condition2)) %>%
      dplyr::filter(run %in% colnames(sig_counts_scaled)) %>%
      dplyr::arrange(condition, host_subject_id)
    
    if (nrow(meta_filtered) == 0) {
      cat("Filtered metadata has 0 rows - Skipping\n")
      next
    }
    
    sig_counts_scaled <- sig_counts_scaled[, meta_filtered$run, drop = FALSE]
    colnames(sig_counts_scaled) <- meta_filtered$host_subject_id
    
    # Annotations
    unique_hosts <- unique(meta_filtered$host_subject_id)
    host_colors <- setNames(colorRampPalette(brewer.pal(9, "Set1"))(length(unique_hosts)), unique_hosts)
    
    unique_conditions <- unique(meta_filtered$condition)
    if (length(unique_conditions) == 2) {
      condition_colors <- setNames(c("#1f78b4", "#33a02c"), unique_conditions)
    } else {
      condition_colors <- setNames(RColorBrewer::brewer.pal(max(3, length(unique_conditions)), "Set2"), unique_conditions)
    }
    
    annotation_col <- ComplexHeatmap::HeatmapAnnotation(
      Condition = meta_filtered$condition,
      Host_ID = ComplexHeatmap::anno_simple(meta_filtered$host_subject_id, col = host_colors),
      col = list(Condition = condition_colors),
      annotation_legend_param = list(Condition = list(title = "Condition")),
      show_annotation_name = TRUE
    )
    
    heatmap_colors <- circlize::colorRamp2(c(-2, 0, 2), c("#2c7bb6", "#ffffbf", "#d7191c"))
    
    heatmap_plot <- ComplexHeatmap::Heatmap(
      sig_counts_scaled,
      name = "Expression",
      col = heatmap_colors,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 8),
      column_names_gp = grid::gpar(fontsize = 8, col = host_colors[colnames(sig_counts_scaled)]),
      top_annotation = annotation_col,
      column_split = factor(meta_filtered$condition, levels = unique(meta_filtered$condition)),
      heatmap_legend_param = list(title = "Expression", legend_direction = "vertical"),
      border = TRUE
    )
    
    # Final output name with dataset
    outname <- paste0(dataset_name, "_", comparison, ".png")
    png_filename <- file.path(output_dir, outname)
    
    png(png_filename, width = 4800, height = 3000, res = 600)
    ComplexHeatmap::draw(heatmap_plot)
    dev.off()
    
    cat("Heatmap saved:", png_filename, "\n")
  }
}



# List DA files
all_files <- list.files("DA_results/da_sig", pattern = "\\.csv$", full.names = TRUE)

# Filter for dataset-specific files
dieguez_files <- grep("dieguez_", all_files, value = TRUE)
ev_files <- grep("ev_", all_files, value = TRUE)

# Generate heatmaps
ensure_dir("DA_results/heatmaps")
generate_all_heatmaps(dieguez_files, dieguez_complete_uniref_counts, metadata_dieguez, "dieguez", "DA_results/heatmaps")
generate_all_heatmaps(ev_files, ev_complete_uniref_counts, metadata_ev, "ev", "DA_results/heatmaps")


# Add on heatmap for the ones with a lot of genes
# generate_all_heatmaps_1 <- function(da_files, counts_matrix, metadata, dataset_name, output_dir) {
#   if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
#   
#   metadata$condition <- metadata$condition %>%
#     trimws() %>%
#     stringr::str_replace_all("\\s{2,}", " ")
#   
#   metadata_conditions <- unique(metadata$condition)
#   
#   clean_condition <- function(cond_string) {
#     cond_string %>%
#       gsub("_", " - ", .) %>%
#       gsub("\\s{2,}", " ", .) %>%
#       stringr::str_trim() %>%
#       stringr::str_replace_all("\\bFl Ar\\b", "Fl-Ar") %>%
#       stringr::str_replace_all("\\bCI\\b", "CI") %>%
#       stringr::str_replace_all("\\bCF\\b", "CF") %>%
#       tools::toTitleCase()
#   }
#   
#   
#   for (da_file in da_files) {
#     file_base <- basename(da_file)
#     comparison <- file_base %>%
#       gsub("^dieguez_|^ev_", "", .) %>%
#       gsub("\\.csv$", "", .)
#     
#     cat("\nProcessing file:", file_base, "-> Comparison:", comparison, "\n")
#     
#     da_data <- readr::read_csv(da_file, show_col_types = FALSE)
#     significant_taxa <- unique(da_data$Taxon)
#     sig_counts <- counts_matrix[rownames(counts_matrix) %in% significant_taxa, , drop = FALSE]
#     
#     if (nrow(sig_counts) == 0) {
#       cat("No significant genes found for", comparison, "- Skipping heatmap\n")
#       next
#     }
#     
#     sig_counts_log <- log2(sig_counts + 1)
#     sig_counts_scaled <- t(scale(t(sig_counts_log)))
#     
#     conditions <- unlist(strsplit(comparison, "_vs_"))
#     conditions <- lapply(conditions, clean_condition)
#     conditions <- unlist(conditions)
#     
#     if (length(conditions) != 2) {
#       cat("Invalid comparison format (expected 2 conditions):", comparison, "\n")
#       next
#     }
#     
#     cat("Extracted Conditions:", conditions[1], "vs", conditions[2], "\n")
#     cat("Available metadata conditions:\n")
#     print(metadata_conditions)
#     
#     match_condition <- function(cond) {
#       exact <- metadata_conditions[metadata_conditions == cond]
#       if (length(exact) > 0) return(exact[1])
#       fuzzy <- metadata_conditions[agrep(cond, metadata_conditions, ignore.case = TRUE, max.distance = 0.2)]
#       if (length(fuzzy) > 0) return(fuzzy[1])
#       return(NA)
#     }
#     
#     condition1 <- match_condition(conditions[1])
#     condition2 <- match_condition(conditions[2])
#     
#     if (is.na(condition1) || is.na(condition2) || condition1 == condition2) {
#       cat("Could not match both conditions or they are identical:", condition1, condition2, "\n")
#       next
#     }
#     
#     cat("Final Matched Conditions:", condition1, "vs", condition2, "\n")
#     
#     meta_filtered <- metadata %>%
#       dplyr::filter(condition %in% c(condition1, condition2)) %>%
#       dplyr::filter(run %in% colnames(sig_counts_scaled)) %>%
#       dplyr::arrange(condition, host_subject_id)
#     
#     if (nrow(meta_filtered) == 0) {
#       cat("Filtered metadata has 0 rows - Skipping\n")
#       next
#     }
#     
#     sig_counts_scaled <- sig_counts_scaled[, meta_filtered$run, drop = FALSE]
#     colnames(sig_counts_scaled) <- meta_filtered$host_subject_id
#     
#     # Annotations
#     unique_hosts <- unique(meta_filtered$host_subject_id)
#     host_colors <- setNames(colorRampPalette(brewer.pal(9, "Set1"))(length(unique_hosts)), unique_hosts)
#     
#     unique_conditions <- unique(meta_filtered$condition)
#     if (length(unique_conditions) == 2) {
#       condition_colors <- setNames(c("#1f78b4", "#33a02c"), unique_conditions)
#     } else {
#       condition_colors <- setNames(RColorBrewer::brewer.pal(max(3, length(unique_conditions)), "Set2"), unique_conditions)
#     }
#     
#     annotation_col <- ComplexHeatmap::HeatmapAnnotation(
#       Condition = meta_filtered$condition,
#       Host_ID = ComplexHeatmap::anno_simple(meta_filtered$host_subject_id, col = host_colors),
#       col = list(Condition = condition_colors),
#       annotation_legend_param = list(Condition = list(title = "Condition")),
#       show_annotation_name = TRUE
#     )
#     
#     heatmap_colors <- circlize::colorRamp2(c(-2, 0, 2), c("#2c7bb6", "#ffffbf", "#d7191c"))
#     
#     heatmap_plot <- ComplexHeatmap::Heatmap(
#       sig_counts_scaled,
#       name = "Expression",
#       col = heatmap_colors,
#       cluster_rows = TRUE,
#       cluster_columns = FALSE,
#       show_row_names = TRUE,
#       show_column_names = TRUE,
#       row_names_gp = grid::gpar(fontsize = 6),
#       column_names_gp = grid::gpar(fontsize = 8, col = host_colors[colnames(sig_counts_scaled)]),
#       top_annotation = annotation_col,
#       column_split = factor(meta_filtered$condition, levels = unique(meta_filtered$condition)),
#       heatmap_legend_param = list(title = "Expression", legend_direction = "vertical"),
#       border = TRUE
#     )
#     
#     # Final output name with dataset
#     outname <- paste0(dataset_name, "_", comparison, ".png")
#     png_filename <- file.path(output_dir, outname)
#     
#     png(png_filename, width = 4800, height = 8000, res = 600)
#     ComplexHeatmap::draw(heatmap_plot)
#     dev.off()
#     
#     cat("Heatmap saved:", png_filename, "\n")
#   }
# }
# 
# 
# 
# # List DA files
# all_files <- list.files("DA_results/da_sig", pattern = "\\.csv$", full.names = TRUE)
# 
# # Filter for dataset-specific files
# dieguez_files <- grep("dieguez_", all_files, value = TRUE)
# ev_files <- grep("ev_", all_files, value = TRUE)
# 
# # Generate heatmaps
# ensure_dir("DA_results_1/heatmaps")
# generate_all_heatmaps_1(ev_files, ev_complete_uniref_counts, metadata_ev, "ev", "DA_results_1/heatmaps")
