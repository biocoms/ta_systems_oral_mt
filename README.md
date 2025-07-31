# Toxin–Antitoxin Systems in the Oral Microbiome

This repository accompanies the computational and analytical framework used to characterize **toxin–antitoxin (TA) systems** in the oral microbiome using metatranscriptomic data from two publicly available cohorts: **Dieguez et al.** and **Ev** datasets.

The analysis integrates differential expression, curated functional annotations, and detailed visualizations to identify and interpret the transcriptional activity of TA gene pairs in healthy and caries-associated and caries-treated oral microbiomes. All scripts, intermediate data files, and processed outputs are included to enable reproducibility, transparency, and downstream reuse.

${\color{red}Disclaimer}$ - The raw FASTQ reads, intermediate preprocessing files, and full functional annotations are not included here due to storage constraints. However, detailed scripts and a step-by-step tutorial are provided below.

---

- [Toxin–Antitoxin Systems in the Oral Microbiome](#toxinantitoxin-systems-in-the-oral-microbiome)
  - [Repository Overview \& Structure](#repository-overview--structure)
  - [Installation instructions](#installation-instructions)
    - [Clone the Repository](#clone-the-repository)
    - [Install packages and download databases](#install-packages-and-download-databases)
  - [Metatranscriptomic Data processing Pipeline](#metatranscriptomic-data-processing-pipeline)
    - [Download Raw Reads (SRA/ENA)](#download-raw-reads-sraena)
    - [Quality Control and Trimming (Trim Galore)](#quality-control-and-trimming-trim-galore)
    - [Host and rRNA Removal (SortMeRNA)](#host-and-rrna-removal-sortmerna)
  - [Functional Profiling (HUMAnN)](#functional-profiling-humann)
    - [Directory Structure](#directory-structure)
  - [Functional annotations](#functional-annotations)
  - [UniRef90 ID Extraction and FASTA Retrieval](#uniref90-id-extraction-and-fasta-retrieval)
    - [Database annotations](#database-annotations)
  - [Downstream Analysis: Differential Expression and Visualization](#downstream-analysis-differential-expression-and-visualization)
    - [Environment Setup](#environment-setup)
    - [Data Cleaning and Preprocessing](#data-cleaning-and-preprocessing)
    - [UniRef90 Toxin Cluster Overlap and Venn Analysis](#uniref90-toxin-cluster-overlap-and-venn-analysis)
    - [Differential Abundance Analysis (ANCOM-BC)](#differential-abundance-analysis-ancom-bc)
    - [Volcano Plot Visualization](#volcano-plot-visualization)
    - [Expression Heatmaps (ComplexHeatmap)](#expression-heatmaps-complexheatmap)
    - [Gene Set Overlap Analysis: UpSet and Venn Panels](#gene-set-overlap-analysis-upset-and-venn-panels)
      - [UpSet Plot of Significant TA Genes](#upset-plot-of-significant-ta-genes)
      - [Comparison Label Harmonization](#comparison-label-harmonization)
      - [Venn Diagram Panels (Three-Way Intersection)](#venn-diagram-panels-three-way-intersection)
      - [Venn Panel Reference Table](#venn-panel-reference-table)
  - [Functional Annotations](#functional-annotations-1)
    - [Parsing TADB3 DIAMOND Results](#parsing-tadb3-diamond-results)
    - [Annotating Significant Genes with TADB3](#annotating-significant-genes-with-tadb3)
    - [Annotating Significant Genes with VFDB](#annotating-significant-genes-with-vfdb)
    - [Annotating Significant Genes with eggNOG-mapper](#annotating-significant-genes-with-eggnog-mapper)
    - [Annotating Significant Genes with InterProScan](#annotating-significant-genes-with-interproscan)
  - [Validated TA Gene Clusters](#validated-ta-gene-clusters)
    - [Summary Statistics by Condition](#summary-statistics-by-condition)
    - [Condition-Wise Expression Barplots](#condition-wise-expression-barplots)
  - [Authors and Maintainers](#authors-and-maintainers)
  - [Questions and Feedback](#questions-and-feedback)

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

${\color{blue}Note}$ - All files larger than 50MB are tracked using Git LFS. Access instructions are provided below.

---

## Installation instructions

All analyses were developed and executed on a Linux-based HPC cluster and assume working familiarity with common bioinformatics tools plus basic R/Python scripting.

### Clone the Repository

```bash

git clone https://github.com/biocoms/ta_systems_oral_mt.git
cd ta_systems_oral_mt

```

### Install packages and download databases

| Tool              | Version        | Source / Install Method                                                     | Documentation Link                                                                  |
| ----------------- | -------------- | --------------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| **FastQC**        | `0.12.1`       | `conda -c bioconda`                                                         | [FastQC Docs](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)           |
| **Cutadapt**      | `2.10`         | `conda -c bioconda`                                                         | [Cutadapt Docs](https://cutadapt.readthedocs.io/en/stable/)                         |
| **Trim Galore**   | `0.6.10`       | `conda -c bioconda`                                                         | [Trim Galore Docs](https://github.com/FelixKrueger/TrimGalore)                      |
| **SortMeRNA**     | `4.3.7`        | `conda -c bioconda`                                                         | [SortMeRNA GitHub](https://github.com/biocore/sortmerna)                            |
| **HUMAnN**        | `3.9`          | `conda -c bioconda`                                                         | [HUMAnN Docs](https://github.com/biobakery/humann)                          |
| **Bowtie2**       | `>=2.4` (auto) | `conda -c bioconda`                                                         | [Bowtie2 Docs](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)              |
| **DIAMOND**       | `>=2.1` (auto) | `conda -c bioconda`                                                         | [DIAMOND Docs](https://github.com/bbuchfink/diamond)                                |
| **eggNOG-mapper** | `2.1.12`       | `conda -c bioconda`                                                         | [eggNOG-mapper Docs](https://github.com/eggnogdb/eggnog-mapper)                     |
| **HMMER**         | `auto`         | `conda -c bioconda`                                                         | [HMMER Docs](http://hmmer.org/)                                                     |
| **OpenJDK**       | `11`           | `conda -c conda-forge`                                                      | [OpenJDK Docs](https://openjdk.org/)                                                |
| **InterProScan**  | `5.72-103.0`   | [Manual FTP](https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.72-103.0/) | [InterProScan Docs](https://interproscan-docs.readthedocs.io/)                      |

To ease the hassle of the installation, we have created a bash script for one step installation.
${\color{blue}Note}$ - Ensure this github is cloned, you are inside this folder and the [conda](https://www.anaconda.com/docs/getting-started/miniconda/install) is installed.

```bash
bash scripts/setup.sh
```

If any large files are missing or stubbed out, ensure [Git LFS](https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage) is installed. You can also install `git lfs` via `conda` and activate in `git` as shown below.

```bash

git lfs install
git lfs pull
git lfs track "filename or file extension"
git add .gitattributes

```

## Metatranscriptomic Data processing Pipeline

The following section outlines the complete preprocessing pipeline to go from raw reads to functional gene and pathway profiles using publicly available tools.

### Download Raw Reads (SRA/ENA)

Enter these BioProject numbers in [SRA-Explorer](https://sra-explorer.info), copy the FTP links of the paired-end FASTQ files to a .txt file and download using `wget`

- Dieguez: BioProject [PRJNA712952](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA712952)

- Ev: BioProject [PRJNA930965](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA930965)

```bash

wget -i dieguez.txt
wget -i ev.txt

```

Make sure dieguez.txt and ev.txt each contain full FTP links for both R1 and R2 reads e.g.:

```bash

ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR233/020/SRR23351020/SRR23351020_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR233/020/SRR23351020/SRR23351020_2.fastq.gz
.......
.......

```

### Quality Control and Trimming (Trim Galore)

[Trim Galore](https://github.com/FelixKrueger/TrimGalore) is a wrapper tool that combines FastQC and Cutadapt to perform quality filtering and adapter trimming.

Outputs:

- Trimmed paired FASTQ files (e.g., _1_trimmed.fastq.gz,_2_trimmed.fastq.gz)
- FastQC quality reports (_fastqc.html,_fastqc.zip)

Command Example:

```bash

trim_galore --paired --fastqc --gzip --phred33 --length 50 \
  --output_dir dieguez/trimmed_reads sample_1.fastq.gz sample_2.fastq.gz

```

- --paired: Indicates paired-end reads.

- --fastqc: Runs FastQC before and after trimming for QC reports.

- --gzip: Compresses the output files in .gz format.

- --phred33: Specifies base quality encoding (standard for Illumina).
  
- --length 50: Discards reads shorter than 50 bp after trimming.
  
- --output_dir: Output folder for trimmed reads and FastQC reports.

This step is performed for (all sequences) both Dieguez and Ev datasets.
FastQC reports will help assess sequence quality pre- and post-trimming.

### Host and rRNA Removal (SortMeRNA)

[SortMeRNA](https://github.com/sortmerna/sortmerna) filters out unwanted rRNA and host (human) reads using curated databases.

Databases Used:

- SILVA rRNA databases (16S, 23S, 5S, etc.)
- Rfam (non-coding RNAs)
- GRCh38 (human genome and transcriptome)

All reference FASTA files are stored in: dbs/sortmerna_databases/

Inputs:

Trimmed paired FASTQ files (trimmed_reads/*_1_trimmed.fastq.gz and trimmed_reads/*_2_trimmed.fastq.gz)

Outputs:

- *_aligned.fastq.gz — reads matching rRNA/host
- *_unaligned.fastq.gz — filtered reads for downstream analysis

Command Template:

```bash

sortmerna --ref dbs/ref_sortmerna/silva-bac-16s-id90.fasta \
          --reads sample_1_trimmed.fastq.gz \
          --reads sample_2_trimmed.fastq.gz \
          --fastx \
          --aligned output_dir/sample_aligned \
          --other output_dir/sample_unaligned \
          --threads 32                        #flexible as per available cores

```

- --ref: Specifies each reference database for filtering.
  
- --reads: Input FASTQ files (both forward and reverse).
  
- --fastx: Ensures FASTQ format is preserved for input/output.
  
- --aligned: Output prefix for reads that matched references.
  
- --other: Output prefix for reads that did not match references.
  
- --threads: Number of threads to use for parallel processing.

Filtered reads are saved into:

- dieguez/trimmed_reads/sortmerna_unaligned/
- ev/trimmed_reads/sortmerna_unaligned/

These unaligned reads are used as input for HUMAnN.

## Functional Profiling (HUMAnN)

[HUMAnN3](https://github.com/biobakery/humann) (The HMP Unified Metabolic Analysis Network) is used to profile gene families and metabolic pathways from host-filtered microbial reads.

Databases Required:

- ChocoPhlAn: nucleotide-level mapping
- UniRef90 or UniRef50: translated protein database for function

Inputs:

- .fq.gz files from sortmerna_unaligned/ directory

Outputs (for each sample):

- {sample}_genefamilies.tsv
- {sample}_pathabundance.tsv
- {sample}_pathcoverage.tsv

Command Template:

```bash

humann -i sample_unaligned.fq.gz \
       -o humann_output/sample_name \
       --threads 32 \
       --nucleotide-database humann/dbs/chocophlan \
       --protein-database humann/dbs/uniref \
       --verbose

```

- -i: Input FASTQ file (filtered reads).
- -o: Output directory for results.
- --threads: Number of parallel threads for speed.
- --nucleotide-database: Path to ChocoPhlAn nucleotide database.
- --protein-database: Path to UniRef protein database.
- --verbose: Prints detailed progress during execution.

Run separately for each dataset (Dieguez and Ev).

Run separately for each dataset (Dieguez and Ev).
This step generates functional profiles that can be used for downstream gene analysis and comparison towards TA gene clusters.

### Directory Structure

```bash
ta_systems_oral_mt/
├── dieguez/
│   ├── trimmed_reads/
│   │   ├── *_1_trimmed.fastq.gz
│   │   ├── *_2_trimmed.fastq.gz
│   │   ├── *_fastqc.html
│   │   ├── *_fastqc.zip
│   │   ├── sortmerna/
│   │   │   ├── <sample>/
│   │   │   │   ├── <sample>_aligned.fastq.gz
│   │   │   │   ├── <sample>_unaligned.fastq.gz
│   │   └── sortmerna_unaligned/
│   │       └── *.fq.gz                     # Filtered (host & rRNA removed) reads
│   └── humann/
│       └── <sample>/
│           ├── <sample>_genefamilies.tsv
│           ├── <sample>_pathabundance.tsv
│           └── <sample>_pathcoverage.tsv
│
├── ev/
│   ├── trimmed_reads/
│   │   ├── *_1_trimmed.fastq.gz
│   │   ├── *_2_trimmed.fastq.gz
│   │   ├── *_fastqc.html
│   │   ├── *_fastqc.zip
│   │   ├── sortmerna/
│   │   │   ├── <sample>/
│   │   │   │   ├── <sample>_aligned.fastq.gz
│   │   │   │   ├── <sample>_unaligned.fastq.gz
│   │   └── sortmerna_unaligned/
│   │       └── *.fq.gz                     # Filtered reads
│   └── humann/
│       └── <sample>/
│           ├── <sample>_genefamilies.tsv
│           ├── <sample>_pathabundance.tsv
│           └── <sample>_pathcoverage.tsv
│
├── dbs/
│   └── sortmerna_databases/                     # SortMeRNA reference databases
│       ├── *.fasta
│   └── humann_databases/
│       ├── chocophlan/                    # HUMAnN nucleotide database
│       └── uniref/                        # HUMAnN protein database
│
├── dieguez.txt                            # FTP links for Dieguez dataset
├── ev.txt                                 # FTP links for Ev dataset
├── ta_process.sh                          # Main preprocessing script
└── log/
    └── processing_pipeline.log            # Cluster job log output

```

---

## Functional annotations

## UniRef90 ID Extraction and FASTA Retrieval

We first extracted unique **UniRef90 gene family IDs** from the `*_genefamilies.tsv` outputs of HUMAnN. These IDs were then used to fetch FASTA sequences for downstream functional annotation using InterProScan, eggNOG-mapper, TADB3, and VFDB.

Inputs:

- Directory of `*_genefamilies.tsv` files from HUMAnN for each sample  
  (`raw_data/dieguez_genefamilies.tsv`, `raw_data/ev_genefamilies.tsv`)

Steps

1. **Extract UniRef90 IDs**
   - Parsed each genefamilies file to collect only those gene families starting with `UniRef90_`.
   - Saved as one ID list per sample in:  
     `uniref_ids/extracted_ids/*_uniref90_ids.txt`
2. **Download FASTA Sequences**
   - Queried the UniProt REST API in chunks (default = 500 IDs/query).
   - Fetched `.fasta` sequences corresponding to UniRef90 IDs per sample.
   - Logs were saved for each sample to monitor download status.

Script:

- `scripts/uniref_idmapping.py`

```python

python scripts/uniref_idmapping.py \
  --input_dir raw_data/dieguez_genefamilies.tsv \
  --output_dir dieguez_uniref_mapped/

```

```python

python scripts/uniref_idmapping.py \
  --input_dir raw_data/ev_genefamilies.tsv \
  --output_dir ev_uniref_mapped/

```

Outputs:

- `uniref_ids/extracted_ids/*.txt` – Sample-wise UniRef90 ID lists
- `uniref_ids/fastas/*.fasta` – Fetched amino acid sequences
- `uniref_ids/logs/*.log` – Chunk-wise log of downloads per sample

### Database annotations

The shell script `database_annotations.sh` automates batch functional annotation using four major tools: **VFDB**, **TADB3**, **eggNOG-mapper**, and **InterProScan**, applied to UniRef90-mapped FASTA sequences from the Dieguez and EV datasets.

Inputs

- **FASTA Sequences**: Protein-coding genes extracted using UniRef90 IDs from significant TA genes.
  - Located in:
    - `dieguez_uniref_mapped/fastas/`
    - `ev_uniref_mapped/fastas/`

- **Databases**:
  - `vfdb/VFDB_prot.dmnd`
  - `tadb3/tadb3_combined.fasta/*.dmnd`
  - `dbs/eggnog_data/eggnog_proteins.dmnd`
  - `dbs/interproscan-5.72-103.0/interproscan.sh`

- **Threads**: VFDB/TADB3 (64), eggNOG (80), InterProScan uses all available CPUs.

Steps

1. **VFDB (DIAMOND Blastp)**  
   Matches each UniRef90 FASTA file to the VFDB protein database.  
   Output format: tab-delimited DIAMOND blastp (`.tsv`)

2. **TADB3 (DIAMOND Blastp to Multiple DBs)**  
   Iteratively matches each FASTA to all `*.dmnd` toxin/antitoxin databases from TADB3.

3. **eggNOG-mapper**  
   Runs `emapper.py` with:
   - Minimum 50% identity (`--pident`)
   - Minimum 80% query coverage (`--query_cover`)
   - DIAMOND mode, using custom eggNOG directory

4. **InterProScan**  
   Annotates each FASTA using InterProScan with GO term and IPR lookups.  
   Output format: tab-separated `.tsv` with all domain and functional annotations.

Outputs

| Tool           | Dieguez Output Directory           | EV Output Directory                |
|----------------|-------------------------------------|-------------------------------------|
| **VFDB**       | `vfdb/dieguez/`                    | `vfdb/ev/`                          |
| **TADB3**      | `tadb3/dieguez/`            | `tadb3/ev/`                  |
| **eggNOG**     | `eggnog/dieguez/`         | `eggnog/ev/`               |
| **InterProScan** | `interproscan/dieguez/`   | `interproscan/ev/`          |

Execution

To run all annotations:

```bash

bash database_annotations.sh

```

---

## Downstream Analysis: Differential Expression and Visualization

All downstream analyses were conducted in R using the `ta_systems_oral_mt.Rproj` environment with strict reproducibility ensured via `renv`.

### Environment Setup

Activate the project-local R environment using:

```r
# From within the R project root
install.packages("renv")
renv::restore()
```

This installs all required packages at pinned versions as declared in `renv.lock`.

Alternatively, for manual setup, install the following CRAN and Bioconductor packages:

```r

# CRAN packages
install.packages(c(
  "tidyverse", "dplyr", "ggplot2", "gridExtra", "patchwork", 
  "VennDiagram", "RColorBrewer", "pheatmap", "UpSetR"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "phyloseq", "ANCOMBC", "ComplexHeatmap", "circlize", "microbiome"
))

```

### Data Cleaning and Preprocessing

Raw gene family expression files (`*_genefamilies.tsv`) generated by HUMAnN were cleaned and normalized as follows:

- Retained only gene families prefixed with `UniRef90_`
- Removed unmapped or low-abundance entries (e.g., `UNMAPPED`)
- Standardized sample identifiers by removing `.RPKs` suffixes
- Ensured uniqueness of gene family rows

Output Files:

- `processed_data/processed_dieguez_genes.csv`
- `processed_data/processed_ev_genes.csv`

### UniRef90 Toxin Cluster Overlap and Venn Analysis

We intersected processed gene family tables with the curated UniRef90 toxin cluster list (`uniref90_toxin_toxic.tsv`) to identify shared and unique transcriptional signals across datasets.

- Computed overlaps between Dieguez and EV gene families with UniRef90 toxin genes
- Separated results into:
  - Dieguez-only overlaps
  - EV-only overlaps
  - Shared UniRef90-positive genes

Output Matrices:

- `uniref_tox_abundance/dieguez_final_uniref90_mapped.csv`
- `uniref_tox_abundance/ev_final_uniref90_mapped.csv`
- `uniref_tox_abundance/intersection_dieguez_final_uniref90_mapped.csv`
- `uniref_tox_abundance/intersection_ev_final_uniref90_mapped.csv`

Venn Diagram

- `venn_diagrams/venn_dieguez_ev_uniref.png`

### Differential Abundance Analysis (ANCOM-BC)

Differential abundance testing was conducted using [ANCOM-BC](https://github.com/FrederickHuangLin/ANCOMBC) independently within each dataset.

Workflow:

- Constructed `phyloseq` objects from UniRef90-filtered matrices
- Enumerated all valid pairwise condition comparisons
- Applied `ancombc()` with the model formula `~ condition`
- Exported full and filtered (significant) result tables

Output Tables:

- `DA_results/da_complete/*.csv`  
- `DA_results/da_sig/*.csv`

### Volcano Plot Visualization

For each ANCOM-BC comparison, volcano plots were generated to highlight gene-level statistical significance and fold-change magnitude.

- **x-axis:** $\log_{2}$ fold-change (`lfc`)
- **y-axis:** –\log_{10}$(p-value)
- **Color scheme:** red (significant), gray (non-significant)
- Vertical and horizontal thresholds marked at ±1 $\log_{2}$(FC) and p = 0.05

Output Directory:

- `DA_results/volcano_plot/*.png`

### Expression Heatmaps (ComplexHeatmap)

Significantly differentially expressed TA gene clusters were visualized using scaled heatmaps.

- Expression values were $\log_{2}$-transformed with pseudocounts
- Row-wise scaling was applied to highlight relative expression shifts
- Columns were annotated using condition and `host_subject_id`
- Color palettes were customized to reflect biological and sample group distinctions

Output Directory:

- `DA_results/heatmaps/*.png`

### Gene Set Overlap Analysis: UpSet and Venn Panels

To summarize and visualize shared transcriptional activity of UniRef90-mapped toxin–antitoxin genes across conditions, we generated both **UpSet plots** and **Venn diagrams**.

#### UpSet Plot of Significant TA Genes

The UpSet plot captures the intersection of all significantly expressed UniRef90 genes across validated pairwise comparisons.

Highlights:

- **Rows:** Gene clusters
- **Columns:** Comparisons (Dieguez and EV)
- **Bars:** Number of shared genes across comparison sets
- **Queries:**
  - Red: genes annotated as **toxins**
  - Green: genes annotated as **antitoxins**

Color Mapping:

- Dieguez comparisons: `#E82561`
- EV comparisons: `purple`

Output:

- `DA_results/upset_plot/upset_plot.png`

#### Comparison Label Harmonization

To ensure consistency in comparison labels across datasets, condition names were programmatically cleaned and standardized. Reverse contrasts (e.g., `B vs A` and `A vs B`) were collapsed.

#### Venn Diagram Panels (Three-Way Intersection)

Nine thematic panels (B to J) were designed to illustrate 3-way overlap across clinically meaningful comparisons.

Examples:

- **Panel B:** Caries–CA* vs Healthy–CF (AAa, CAi, CAs)
- **Panel H:** Healthy-only timepoint comparisons
- **Panel J:** Healthy vs Caries across all stages

Each Venn diagram shows:

- Size-labeled sets with word-wrapped condition labels
- Color-coded ellipses from `RColorBrewer::Set2`
- Gene counts in intersections and unique sections

Output Directory:

- `DA_results/venn_panels/`

File Naming Convention:

- `Venn_Panel_<Letter>.png`

**Supporting Tables:**

- `All_Venn_Sets_Combined.csv`: All individual sets per panel
- `All_Venn_Intersections_3way.csv`: Genes shared by all three comparisons per panel

#### Venn Panel Reference Table

| Panel | Comparisons |
|-------|-------------|
| B     | Caries-CAa vs Healthy-CF, Caries-CAi vs Healthy-CF, Caries-CAs vs Healthy-CF |
| C     | Caries-CAa vs Healthy-CI, Caries-CAi vs Healthy-CI, Caries-CAs vs Healthy-CI |
| D     | Healthy-CF vs Healthy-CI, Caries-CAa vs Healthy-CF, Caries-CAa vs Healthy-CI |
| E     | Healthy-CF vs Healthy-CI, Caries-CAi vs Healthy-CF, Caries-CAi vs Healthy-CI |
| F     | Healthy-CF vs Healthy-CI, Caries-CAs vs Healthy-CF, Caries-CAs vs Healthy-CI |
| G     | Caries-CAa vs Caries-CAi, Caries-CAa vs Caries-CAs, Caries-CAi vs Caries-CAs |
| H     | Healthy-Baseline vs Healthy-Fl, Healthy-Baseline vs Healthy-Fl-Ar, Healthy-Fl vs Healthy-Fl-Ar |
| I     | Caries-Baseline vs Caries-Fl, Caries-Baseline vs Caries-Fl-Ar, Caries-Fl vs Caries-Fl-Ar |
| J     | Caries-Baseline vs Healthy-Baseline, Caries-Fl vs Healthy-Fl, Caries-Fl-Ar vs Healthy-Fl-Ar |

## Functional Annotations

### Parsing TADB3 DIAMOND Results

`tadb3_master_table.py` parses DIAMOND alignment outputs against the TADB3 database and generates a consolidated annotation table for all matched genes across Dieguez and Ev datasets.

Inputs:

- DIAMOND outputs:
  - `tadb3/diamond_dieguez/*.tsv`
  - `tadb3/diamond_ev/*.tsv`
- Reference FASTA:
  - `tadb3/tadb3_combined.fasta`

Steps:

1. Parse TADB3 FASTA to extract `Protein_Accession`, `Genome_Accession`, and `Organism`.
2. Read DIAMOND result files per sample; each file contains UniRef90 hits to TADB3 proteins.
3. Merge DIAMOND hits with metadata extracted from the FASTA.
4. Annotate each hit with:
   - `Toxin_Antitoxin`: inferred from ID patterns (`_tox_`, `_antitox_`)
   - `Validation_Type`: `Experimental` or `Computational` based on UniRef90 ID
5. Record `Sample_Name` and dataset (`dieguez` or `ev`) for each entry.

Output:

- `tadb3/tadb3_tox_antitox_main_table_1.tsv`: master table with annotation metadata for all TADB3-aligned genes across both datasets.

### Annotating Significant Genes with TADB3

`tadb3_sig_genes.py` filters and annotates TADB3-aligned genes that are differentially expressed in the Dieguez or EV datasets. This script integrates DIAMOND-based hits with metadata from the TADB3 master table and augments entries with protein descriptions retrieved from NCBI.

Inputs:

- `DA_results/sig_genes.txt`: List of significant UniRef90 gene IDs
- `tadb3/master_table_tadb3.txt`: Annotated master table from `tadb3_master_table.py`

Steps:

1. Load significant UniRef90 IDs.
2. Subset TADB3 hits for these IDs and retain only high-confidence matches (Bit Score ≥ 90).
3. Standardize organism names to species-level.
4. Query NCBI Entrez for full protein descriptions using `Protein_Accession` and clean the output.
5. For each gene, perform annotation collapse (e.g., most common `Organism`, validation type, aggregated annotations).

Output:

- `tadb3/sig_genes_annotations.csv`: Final annotation table for significant TA genes enriched with curated metadata and Entrez descriptions.

### Annotating Significant Genes with VFDB

`vfdb_sig_genes.py` filters DIAMOND alignments against the VFDB (Virulence Factor Database) to extract and annotate significant genes identified from the Dieguez and EV datasets.

Inputs:

- `DA_results/sig_genes.txt`: List of significant UniRef90 gene IDs.
- DIAMOND outputs:
  - `vfdb/dieguez/*_vs_VFDB.tsv`
  - `vfdb/ev/*_vs_VFDB.tsv`

Steps:

1. Load significant UniRef90 IDs and subset the DIAMOND VFDB alignments to retain only matches with these genes.
2. Extract GenBank accession numbers from `Query_ID` using regex.
3. Query NCBI Entrez using `Bio.Entrez` to fetch gene names and protein descriptions for each accession.
4. Append sample name, dataset, gene name, and description to each annotated entry.
5. Collapse annotations by gene using custom aggregation logic.

Output:

- `vfdb/sig_genes_annotations.csv`: Annotated table of significant genes aligned to VFDB with NCBI-derived metadata.

### Annotating Significant Genes with eggNOG-mapper

`eggnog_sig_genes.py` extracts and summarizes functional annotations for significant UniRef90 gene families using EggNOG-mapper outputs from the Dieguez and EV datasets.

Inputs:

- `DA_results/sig_genes.txt`: List of significant UniRef90 gene IDs (dot-truncated)
- Annotated `.emapper.annotations` files from EggNOG-mapper:
  - `eggnog/dieguez/*_eggnog.emapper.annotations`
  - `eggnog/ev/*_eggnog.emapper.annotations`

Steps:

1. Parse each `.emapper.annotations` file while ignoring comment lines and correctly identifying the column header.
2. Filter rows where the query matches a significant UniRef90 gene ID.
3. Annotate each record with dataset and sample origin.
4. Collapse annotations per UniRef90 ID by computing the most frequent (mode) entry for each functional field (e.g., `KEGG_ko`, `EC`, `eggNOG_OGs`, `GO_terms`).
5. Track how many and which samples contributed to each gene’s annotation.

Output:

- `eggnog/sig_genes_annotations.csv`: Collapsed annotation table of significant genes, including sample metadata and representative functional terms.

### Annotating Significant Genes with InterProScan

`interproscan_sig_genes.py` parses and annotates InterProScan `.tsv` results for significant UniRef90 gene families, enriching each entry with GO term descriptions and ontology categories.

Inputs:

- `DA_results/sig_genes.txt`: List of significant UniRef90 gene IDs.
- InterProScan output files:
  - `interproscan/dieguez/*.tsv`
  - `interproscan/ev/*.tsv`

Steps:

1. Load `.tsv` files and assign dataset/sample labels.
2. Filter records based on significant gene IDs.
3. Extract GO identifiers (`GO:XXXXXXX`) from the GO column.
4. Query the [EBI QuickGO API](https://www.ebi.ac.uk/QuickGO/) to retrieve term names and ontology categories (BP, MF, CC).
5. Annotate each record with:
   - Full GO term descriptions
   - Separated Biological Process, Molecular Function, Cellular Component
6. Aggregate and collapse annotations for each gene across samples and datasets.

Output:

- `interproscan/sig_genes_annotations.csv`: Clean annotation table with GO terms and structured categories for significant genes.

## Validated TA Gene Clusters

This section summarizes and visualizes the expression behavior of experimentally validated or computationally predicted toxin–antitoxin (TA) system genes based on curated annotations from TADB3. Two components are included:

### Summary Statistics by Condition

The script `valid_gene_summary.R` compiles condition-wise expression summaries for all validated TA system genes across the Dieguez and Ev datasets.

Inputs:

- `processed_data/processed_dieguez_genes.csv`
- `processed_data/processed_ev_genes.csv`
- `raw_data/metadata_dieguez.csv`
- `raw_data/metadata_ev.csv`
- `tadb3/tadb3_tox_antitox_main_table.tsv`

Steps:

1. Subset expression matrices to retain only TA genes annotated in the TADB3 table.
2. Merge with metadata to assign each sample to a condition group.
3. Apply a pseudocount (0.01) and compute `log10(expression + 0.01)` per sample.
4. For each gene × condition × dataset, compute:
   - `n`, `non-zero sample count`
   - `min`, `max`, `mean`, `median`
   - `mean_log10`, `median_log10`
   - Gene label (`Toxin` or `Antitoxin`)
5. Merge Dieguez and EV summaries into a combined summary table.

Output:

- `DA_results/debug_gene_condition_summary_combined.csv`

### Condition-Wise Expression Barplots

The script `barplots.R` creates annotated barplots for each significant TA gene identified in `sig_genes_combined.csv`, stratified by dataset and gene type.

Inputs:

- `processed_data/processed_dieguez_genes.csv`
- `processed_data/processed_ev_genes.csv`
- `raw_data/metadata_dieguez.csv`
- `raw_data/metadata_ev.csv`
- `DA_results/sig_genes_combined.csv`
- `tadb3/tadb3_tox_antitox_main_table.tsv`

Steps:

1. Subset significant genes by type (`Toxin`, `Antitoxin`).
2. Plot per-gene barplots using:
   - Y-axis: `log1p(expression)`
   - X-axis: condition group (ordered)
   - Mean ± standard error
3. Overlay ANCOM-BC adjusted p-values as:
   - Brackets for significant comparisons only
   - Significance stars in black
4. Customize plot aesthetics:
   - Red text for toxins, green text for antitoxins
   - One dot per condition per gene (no jitter)

Outputs:

- `DA_results/barplots/toxin_barplots/*.png`
- `DA_results/barplots/antitoxin_barplots/*.png`

## Authors and Maintainers

[Shri Vishalini Rajaram](https://github.com/shrivishalinirajaram), [Priyanka Singh](https://github.com/decoder108) and [Erliang Zeng](https://github.com/zerl)

## Questions and Feedback

This repository is maintained by the Bioinformatics and Computational Systems Biology Lab at the University of Iowa. 

If you notice any errors, have suggestions for improvement, or simply want to discuss aspects of the analysis, we genuinely welcome your feedback. Please feel free to open an issue or contact us. We are always happy to engage with the community.
