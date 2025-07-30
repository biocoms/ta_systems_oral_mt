# Toxin–Antitoxin Systems in the Oral Microbiome

This repository accompanies the computational and analytical framework used to characterize **toxin–antitoxin (TA) systems** in the oral microbiome using metatranscriptomic data from two publicly available cohorts: **Dieguez et al.** and **Ev** datasets.

The analysis integrates differential expression, curated functional annotations, and detailed visualizations to identify and interpret the transcriptional activity of TA gene pairs in healthy and caries-associated and caries-treated oral microbiomes. All scripts, intermediate data files, and processed outputs are included to enable reproducibility, transparency, and downstream reuse.

${\color{red}Disclaimer}$ - All the raw fastq reads, their preprocessing, and functional annotation is not provided due to the storage efficiency. However, detailed scripts and step by step tutorial is mentioned below.

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

The project was developed and tested in a Linux-based environment (HPC cluster), and assumes working familiarity with standard bioinformatics tooling and R/Python scripting.

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

Dieguez: BioProject [PRJNA712952](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA712952&o=acc_s%3Aa)
Ev: BioProject [PRJNA930965](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA930965&o=acc_s%3Aa)

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
[Trim Galore](https://github.com/FelixKrueger/TrimGalore) is a wrapper tool that combines Cutadapt and FastQC to perform quality filtering and adapter trimming.

Outputs:
Trimmed paired FASTQ files (e.g., _1_trimmed.fastq.gz)
FastQC quality reports (*_fastqc.html)

Command Example:

```bash

trim_galore --paired --fastqc --gzip --phred33 --length 50 \
  --output_dir Dieguez/trimmed_reads sample_1.fastq.gz sample_2.fastq.gz

```
This step is performed for both Dieguez and Ev datasets.
FastQC reports will help assess sequence quality pre- and post-trimming.


### 3. Host and rRNA Removal (SortMeRNA)
[SortMeRNA](https://github.com/sortmerna/sortmerna) filters out unwanted rRNA and host (human) reads using curated databases.

Databases Used:

SILVA rRNA databases (16S, 23S, 5S, etc.)
Rfam (non-coding RNAs)
GRCh38 (human genome and transcriptome)

All reference FASTA files are stored in: dbs/ref_sortmerna/

Inputs:

Trimmed paired FASTQ files

Outputs:

*_aligned.fastq.gz — reads matching rRNA/host
*_unaligned.fastq.gz — filtered reads for downstream analysis

Command Template:
```bash
sortmerna --ref dbs/ref_sortmerna/silva-bac-16s-id90.fasta \
          --reads sample_1_trimmed.fastq.gz \
          --reads sample_2_trimmed.fastq.gz \
          --fastx \
          --aligned output_dir/sample_aligned \
          --other output_dir/sample_unaligned \
          --threads 32
```
Filtered reads are saved into:
Dieguez/trimmed_reads/sortmerna_unaligned/
Ev/trimmed_reads/sortmerna_unaligned/

## Functional Profiling (HUMAnN)
[HUMAnN](https://github.com/biobakery/humann) (The HMP Unified Metabolic Analysis Network) is used to profile gene families and metabolic pathways from host-filtered microbial reads.

Databases Required:

ChocoPhlAn: nucleotide-level mapping (download from HUMAnN DBs)
UniRef90 or UniRef50: translated protein database for function

Inputs:

.fq.gz files from sortmerna_unaligned/ directory

Outputs (for each sample):

*_genefamilies.tsv
*_pathabundance.tsv
*_pathcoverage.tsv

Command Template:
```bash
humann -i sample_unaligned.fq.gz \
       -o humann_output/sample_name \
       --threads 32 \
       --nucleotide-database humann/dbs/chocophlan \
       --protein-database humann/dbs/uniref \
       --verbose
```
Run separately for each dataset (Dieguez and Ev).

## Directory Structure
```bash
ta_systems_oral_mt/
├── Dieguez/
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
├── Ev/
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
│   └── ref_sortmerna/                     # SortMeRNA reference databases
│       ├── *.fasta
│
├── humann/
│   └── dbs/
│       ├── chocophlan/                    # HUMAnN nucleotide database
│       └── uniref/                        # HUMAnN protein database
│
├── dieguez.txt                            # FTP links for Dieguez dataset
├── ev.txt                                 # FTP links for Ev dataset
├── ta_process.sh                          # Main preprocessing script
└── log/
    └── processing_pipeline.log            # Cluster job log output
```

## Authors and Maintainers

[Shri Vishalini Rajaram](https://github.com/shrivishalinirajaram), [Priyanka Singh](https://github.com/decoder108) and [Erliang Zeng](https://github.com/zerl)

## Questions and Feedback

This repository is maintained by PhD students actively learning and working at the intersection of microbiome science and computational biology. While we strive for accuracy and reproducibility, this work is part of our training, and mistakes are possible.

If you notice any errors, have suggestions for improvement, or simply want to discuss aspects of the analysis, we genuinely welcome your feedback. Please feel free to open an issue or contact us. We are always learning and happy to engage with the community.
