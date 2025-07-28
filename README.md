# Toxin–Antitoxin Systems in the Oral Microbiome

This repository accompanies the computational and analytical framework used to characterize **toxin–antitoxin (TA) systems** in the oral microbiome using metatranscriptomic data from two publicly available cohorts: **Dieguez et al.** and **Ev** datasets.

The analysis integrates differential expression, curated functional annotations, and detailed visualizations to identify and interpret the transcriptional activity of TA gene pairs in healthy and caries-associated and caries-treated oral microbiomes. All scripts, intermediate data files, and processed outputs are included to enable reproducibility, transparency, and downstream reuse.

${\color{red}\bold{Disclaimer}}$ - All the raw fastq reads, their preprocessing, and functional annotation is not provided due to the storage efficiency. However, detailed scripts and step by step tutorial is mentioned below.

---

## Repository Overview & Structure

The directory structure:

```bash

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

${\color{blue} \bold{Note}}$ - All files larger than 50MB are tracked using Git LFS. Access instructions are provided below.

---

## Installation instructions

The project was developed and tested in a Linux-based environment (HPC cluster), and assumes working familiarity with standard bioinformatics tooling and R/Python scripting.

### 1. Clone the Repository

```bash

git clone https://github.com/biocoms/ta_systems_oral_mt.git
cd ta_systems_oral_mt

```

If any large files are missing or stubbed out, ensure [Git LFS](https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage) is installed. You can also install `git lfs` via `conda` and activate in `git` as shown below.

```bash

git lfs install
git lfs pull

```

### 2. Configure Conda Channels (Highly recommended)

Ensure your Conda installation is correctly configured with strict channel priority.

```bash

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

```

### 3. Create Conda Environment

We recommend creating a dedicated Conda environment (e.g., oral_ta_env) for running preprocessing and database annotation steps.

```bash

conda create -n ta_env python=3.10 -y
conda activate ta_env

```

## Authors and Maintainers

[Shri Vishalini Rajaram](https://github.com/shrivishalinirajaram) and [Priyanka Singh](https://github.com/decoder108)

## Questions and Feedback

This repository is maintained by PhD students actively learning and working at the intersection of microbiome science and computational biology. While we strive for accuracy and reproducibility, this work is part of our training, and mistakes are possible.

If you notice any errors, have suggestions for improvement, or simply want to discuss aspects of the analysis, we genuinely welcome your feedback. Please feel free to open an issue or contact us. We are always learning and happy to engage with the community.
