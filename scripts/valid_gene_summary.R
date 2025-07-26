library(tidyverse)
library(readr)

# Load Data 
counts_dieguez <- read.csv("processed_data/processed_dieguez_genes.csv", row.names = 1, check.names = FALSE)
metadata_dieguez <- read_csv("raw_data/metadata_dieguez.csv") %>%
  select(run, condition) %>%
  mutate(
    run = as.character(run),
    condition = str_replace_all(str_trim(condition), " ?- ?", "-")
  )

counts_ev <- read.csv("processed_data/processed_ev_genes.csv", row.names = 1, check.names = FALSE)
metadata_ev <- read_csv("raw_data/metadata_ev.csv") %>%
  select(run, condition) %>%
  mutate(
    run = as.character(run),
    condition = str_replace_all(str_trim(condition), " ?- ?", "-")
  )

# Load TADB3 and Define Gene Sets 
tadb3 <- read_tsv("tadb3/tadb3_tox_antitox_main_table.tsv")

toxins <- unique(tadb3 %>% filter(Toxin_Antitoxin == "Toxin") %>% pull(UniRef90_ID))
antitoxins <- unique(tadb3 %>% filter(Toxin_Antitoxin == "Antitoxin") %>% pull(UniRef90_ID))
all_genes <- union(toxins, antitoxins)

# Define Function for Summary 
get_gene_condition_summary <- function(counts_df, metadata_df, dataset_name) {
  pseudo <- 0.01
  counts_df$gene <- rownames(counts_df)
  
  long_df <- counts_df %>%
    filter(gene %in% all_genes) %>%
    pivot_longer(-gene, names_to = "SampleID", values_to = "Raw_Expression") %>%
    left_join(metadata_df, by = c("SampleID" = "run")) %>%
    mutate(
      Log_Expression = log10(Raw_Expression + pseudo),
      Type = case_when(
        gene %in% toxins ~ "Toxin",
        gene %in% antitoxins ~ "Antitoxin",
        TRUE ~ "Unknown"
      )
    )
  
  long_df %>%
    group_by(gene, Type, condition) %>%
    summarise(
      N = n(),
      NonZero = sum(Raw_Expression > 0),
      Min = min(Raw_Expression),
      Max = max(Raw_Expression),
      Mean = mean(Raw_Expression),
      Median = median(Raw_Expression),
      Mean_log10 = mean(Log_Expression),
      Median_log10 = median(Log_Expression),
      Dataset = dataset_name,
      .groups = "drop"
    )
}

# Generate Summaries 
summary_dieguez <- get_gene_condition_summary(counts_dieguez, metadata_dieguez, "Dieguez")
summary_ev <- get_gene_condition_summary(counts_ev, metadata_ev, "Ev")

# Combine and Export 
summary_combined <- bind_rows(summary_dieguez, summary_ev)
write_csv(summary_combined, "DA_results/debug_gene_condition_summary_combined.csv")

