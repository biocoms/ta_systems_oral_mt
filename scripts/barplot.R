library(tidyverse)
library(ggtext)
library(ggplot2)
library(patchwork)

# Load Data 
counts_dieguez <- read.csv("processed_data/processed_dieguez_genes.csv", row.names = 1, check.names = FALSE)
counts_ev <- read.csv("processed_data/processed_ev_genes.csv", row.names = 1, check.names = FALSE)
metadata_dieguez <- read_csv("raw_data/metadata_dieguez.csv") %>%
  select(run, condition) %>%
  mutate(condition = str_replace_all(str_trim(condition), " ?- ?", "-"))
metadata_ev <- read_csv("raw_data/metadata_ev.csv") %>%
  select(run, condition) %>%
  mutate(
    condition = str_replace_all(str_trim(condition), " ?- ?", "-"))
tadb3 <- read_tsv("tadb3/tadb3_tox_antitox_main_table.tsv")

# Define Gene Sets 
toxins <- unique(tadb3 %>% filter(Toxin_Antitoxin == "Toxin") %>% pull(UniRef90_ID))
antitoxins <- unique(tadb3 %>% filter(Toxin_Antitoxin == "Antitoxin") %>% pull(UniRef90_ID))

# Fixed Condition Order 
condition_order_dieguez <- c("Healthy-Baseline", "Healthy-Fl", "Healthy-Fl-Ar", "Caries-Baseline", "Caries-Fl", "Caries-Fl-Ar")
condition_order_ev <- c("Healthy-CF", "Healthy-CI", "Caries-CAa", "Caries-CAi", "Caries-CAs")

sig_dieguez <- sig_combined %>%
  filter(str_detect(combinations, "dieguez")) %>%
  mutate(
    comparison_raw = str_remove(combinations, ".*complete_uniref_"),
    group1 = str_replace_all(str_trim(str_split_fixed(comparison_raw, "vs", 2)[, 1]), " ?- ?", "-"),
    group2 = str_replace_all(str_trim(str_split_fixed(comparison_raw, "vs", 2)[, 2]), " ?- ?", "-"),
    gene = Taxon,
    p_val_adj = p.adjust(pval, method = "BH"),
    label = case_when(
      pval < 0.0001 ~ "****",
      pval < 0.001  ~ "***",
      pval < 0.01   ~ "**",
      pval < 0.05   ~ "*",
      TRUE               ~ NA_character_
    )
  ) %>%
  filter(!is.na(label))


sig_ev <- sig_combined %>%
  filter(str_detect(combinations, "^ev_|/ev_|_ev_")) %>%
  mutate(
    comparison_raw = str_remove(combinations, ".*complete_uniref_"),
    group1 = str_replace_all(str_trim(str_split_fixed(comparison_raw, "vs", 2)[, 1]), " ?- ?", "-"),
    group2 = str_replace_all(str_trim(str_split_fixed(comparison_raw, "vs", 2)[, 2])," ?- ?", "-"),
    gene = Taxon,
    p_val_adj = p.adjust(pval, method = "BH"),
    label = case_when(
      pval < 0.0001 ~ "****",
      pval < 0.001  ~ "***",
      pval < 0.01   ~ "**",
      pval < 0.05   ~ "*",
      TRUE               ~ NA_character_
    )
  ) %>%
  filter(!is.na(label))


plot_expression_bar_with_signif <- function(counts, metadata, genes, condition_order, gene_type, sig_df) {
  counts$gene <- rownames(counts)
  pseudo <- 0.01
  
  # Long-form Expression Table 
  df <- counts %>%
    filter(gene %in% genes) %>%
    pivot_longer(-gene, names_to = "SampleID", values_to = "expression") %>%
    left_join(metadata, by = c("SampleID" = "run")) %>%
    filter(!is.na(condition)) %>%
    mutate(
      condition = str_replace_all(condition, " ?- ?", "-"),
      condition = factor(condition, levels = condition_order),
      gene_label = case_when(
        gene_type == "Toxin" ~ paste0("<span style='color:red;'>", gene, "</span>"),
        gene_type == "Antitoxin" ~ paste0("<span style='color:forestgreen;'>", gene, "</span>")
      ),
      log_expr = log1p(expression + pseudo)
    )
  
  # Summary Table (Mean ± SE) 
  df_summary <- df %>%
    group_by(gene, gene_label, condition) %>%
    summarise(
      mean_expr = mean(log_expr, na.rm = TRUE),
      se_expr = sd(log_expr, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Prepare Significance Brackets 
  sig_df_sub <- sig_df %>%
    filter(gene %in% genes) %>%
    mutate(
      gene_label = case_when(
        gene_type == "Toxin" ~ paste0("<span style='color:red;'>", gene, "</span>"),
        gene_type == "Antitoxin" ~ paste0("<span style='color:forestgreen;'>", gene, "</span>")
      )
    ) %>%
    left_join(
      df_summary %>%
        group_by(gene_label) %>%
        summarise(max_expr = max(mean_expr + se_expr, na.rm = TRUE), .groups = "drop"),
      by = "gene_label"
    ) %>%
    group_by(gene_label) %>%
    mutate(y.position = max_expr + 0.3 + 0.3 * row_number()) %>%
    ungroup()
  
  # Plot 
  ggplot(df_summary, aes(x = condition, y = mean_expr, fill = condition)) +
    geom_col(color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = mean_expr - se_expr, ymax = mean_expr + se_expr), width = 0.2) +
    stat_pvalue_manual(
      sig_df_sub,
      label = "label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.3,
      text.size = 3,
      step.increase = 0
    ) +
    facet_wrap(~ gene_label, ncol = 3, scales = "free_y") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Condition", y = "log1p(Mean Expression + 0.01)") +
    theme_bw(base_size = 13) +
    theme(
      strip.text = ggtext::element_markdown(size = 9, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

# Ev toxin/antitoxin bar plots with brackets
p_ev_tox_bar_annot <- plot_expression_bar_with_signif(counts_ev, metadata_ev, toxins, condition_order_ev, "Toxin", sig_ev)
p_ev_ant_bar_annot <- plot_expression_bar_with_signif(counts_ev, metadata_ev, antitoxins, condition_order_ev, "Antitoxin", sig_ev)
p_dieguez_tox_bar_annot <- plot_expression_bar_with_signif(counts_dieguez, metadata_dieguez, toxins, condition_order_dieguez, "Toxin", sig_dieguez)
p_dieguez_ant_bar_annot <- plot_expression_bar_with_signif(counts_dieguez, metadata_dieguez, antitoxins, condition_order_dieguez, "Antitoxin", sig_dieguez)

sig_genes_dieguez <- unique(sig_dieguez$gene)
sig_genes_ev <- unique(sig_ev$gene)

both_sig <- intersect(sig_genes_dieguez, sig_genes_ev)
only_dieguez <- setdiff(sig_genes_dieguez, sig_genes_ev)
only_ev <- setdiff(sig_genes_ev, sig_genes_dieguez)

all_sig_genes <- union(sig_genes_dieguez, sig_genes_ev)


plot_gene_across_datasets <- function(gene, gene_type) {
  plots <- list()
  
  if (gene %in% sig_genes_dieguez) {
    plots[["Dieguez"]] <- plot_expression_bar_with_signif(
      counts_dieguez,
      metadata_dieguez,
      genes = gene,
      condition_order = condition_order_dieguez,
      gene_type = gene_type,
      sig_df = sig_dieguez
    ) +
      ggtitle(paste0(gene, " (Dieguez)"))
  }
  
  if (gene %in% sig_genes_ev) {
    plots[["Ev"]] <- plot_expression_bar_with_signif(
      counts_ev,
      metadata_ev,
      genes = gene,
      condition_order = condition_order_ev,
      gene_type = gene_type,
      sig_df = sig_ev
    ) +
      ggtitle(paste0(gene, " (Ev)"))
  }
  
  # Combine plots (side-by-side or single)
  if (length(plots) == 2) {
    return(plots[[1]] + plots[[2]])
  } else {
    return(plots[[1]])
  }
}
dir.create("DA_results/boxplots/toxin_barplots", showWarnings = FALSE)

for (gene in unique(toxins)) {
  if (gene %in% all_sig_genes) {
    is_in_dieguez <- gene %in% sig_genes_dieguez
    is_in_ev <- gene %in% sig_genes_ev
    num_panels <- as.integer(is_in_dieguez) + as.integer(is_in_ev)
    
    p <- plot_gene_across_datasets(gene, gene_type = "Toxin")
    
    ggsave(
      filename = paste0("DA_results/boxplots/toxin_barplots/", gene, "_barplot.png"),
      plot = p,
      width = 5 * num_panels,
      height = 5,
      dpi = 600
    )
  }
}

dir.create("DA_results/boxplots/antitoxin_barplots", showWarnings = FALSE)

for (gene in unique(antitoxins)) {
  if (gene %in% all_sig_genes) {
    is_in_dieguez <- gene %in% sig_genes_dieguez
    is_in_ev <- gene %in% sig_genes_ev
    num_panels <- as.integer(is_in_dieguez) + as.integer(is_in_ev)
    
    p <- plot_gene_across_datasets(gene, gene_type = "Antitoxin")
    
    ggsave(
      filename = paste0("DA_results/boxplots/antitoxin_barplots/", gene, "_barplot.png"),
      plot = p,
      width = 5 * num_panels,
      height = 5,
      dpi = 600
    )
  }
}

