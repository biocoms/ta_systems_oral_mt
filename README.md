# Toxin–Antitoxin Systems in the Oral Microbiome

This repository accompanies the computational and analytical framework used to characterize **toxin–antitoxin (TA) systems** in the oral microbiome using metatranscriptomic data from two publicly available cohorts: **Dieguez et al.** and **Ev** datasets.

The analysis integrates differential expression, curated functional annotations, and detailed visualizations to identify and interpret the transcriptional activity of TA gene pairs in healthy and caries-associated and caries-treated oral microbiomes. All scripts, intermediate data files, and processed outputs are included to enable reproducibility, transparency, and downstream reuse.

${\color{red} \bold {Disclaimer}}$

All the raw fastq reads, their preprocessing, and functional annotation is not provided due to the storage efficiency. However, detailed scripts and step by step tutorial is mentioned below.

---

## Repository Overview & Structure

The directory structure:

```
ta_systems_oral_mt/
├── DA_results/              # Differential expression results - tables and visualizations
├── eggnog/                  # Functional annotations from eggNOG-mapper
├── interproscan/            # InterProScan annotations for significant genes
├── processed_data/          # Preprocessed and cleaned expression matrices
├── raw_data/                # Raw input files (metadata, UniRef90 toxin-antitoxin files, functional gene tables)
├── scripts/                 # All analysis scripts (R/Python)
├── tadb3/                   # TADB3 toxin–antitoxin annotation tables and sequences
├── uniref_tox_abundance/    # Abundance tables intersecting Dieguez, ev and UniRef90 toxin–antitoxin gene clusters
├── venn_diagrams/           # Venn diagrams used in the manuscript
├── vfdb/                    # VFDB annotations (dieguez/ev-specific)
├── README.md                # Project description and documentation
├── LICENSE                  # License for usage
└── ta_systems_oral_mt.Rproj # R project file for downstream analysis

```
