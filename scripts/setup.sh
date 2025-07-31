#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Configuration
ENV_NAME="ta_env"
INSTALL_DIR="dbs"
SORTMERNA_DB_DIR="$INSTALL_DIR/sortmerna_databases"
INTERPROSCAN_VERSION="5.72-103.0"
INTERPROSCAN_DIR="$INSTALL_DIR/interproscan-${INTERPROSCAN_VERSION}"
INTERPROSCAN_URL="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${INTERPROSCAN_VERSION}"
HUMANN_DB_DIR="$INSTALL_DIR/humann_databases"
EGGNOG_DIR="$INSTALL_DIR/eggnog_data"
LOG_FILE="setup_log.txt"


mkdir -p "$INSTALL_DIR"
touch "$LOG_FILE"
echo "Log file: $LOG_FILE"

echo "Create Conda Environment with Locked Versions" | tee -a "$LOG_FILE"

# Step 1: Configure Channels + Solver
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible  # more compatible
conda config --set auto_update_conda false
conda config --set solver classic  # safer fallback

# Create Environment with Mamba
conda install -n base -c conda-forge mamba -y
mamba create -y -n $ENV_NAME python=3.11

# Activate
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $ENV_NAME

# Install Bioinformatics Tools via Mamba
mamba install -y -c bioconda \
  fastqc=0.12.1 \
  cutadapt=5.0 \
  trim-galore=0.6.10 \
  sortmerna=4.3.7 \
  bowtie2=2.5.4 \
  diamond=2.1.10 \
  mmseqs2 \
  hmmer=3.4 \
  prodigal=2.6.3 \
  eggnog-mapper=2.1.12

mamba install -y -c conda-forge \
  openjdk=11.0.25 \
  git-lfs \
  wget

# Python Package Installs
pip install \
  biopython \
  pandas \
  tqdm \
  requests \
  lxml \
  openpyxl \
  matplotlib \
  seaborn \
  scikit-learn \
  scipy \
  statsmodels
pip install multiqc==1.26
pip install humann==3.9

# Export the env
conda env export --name $ENV_NAME > ta_test.yml

echo "Download and Extract SortMeRNA rRNA Databases" | tee -a "$LOG_FILE"
mkdir -p "$SORTMERNA_DB_DIR"

files=(
  "rfam-5s-database-id98.fasta"
  "silva-bac-16s-id90.fasta"
  "silva-bac-23s-id98.fasta"
  "silva-arc-16s-id95.fasta"
  "silva-arc-23s-id98.fasta"
  "silva-euk-18s-id95.fasta"
  "silva-euk-28s-id98.fasta"
)

BASE_URL="https://raw.githubusercontent.com/sortmerna/sortmerna/master/data/rRNA_databases"

for file in "${files[@]}"; do
  echo "Downloading $file" | tee -a "$LOG_FILE"
  wget -q "$BASE_URL/$file" -O "$SORTMERNA_DB_DIR/$file"
done

echo "Download HUMAnN Databases" | tee -a "$LOG_FILE"
mkdir -p "$HUMANN_DB_DIR"
humann_databases --download chocophlan full "$HUMANN_DB_DIR"
humann_databases --download uniref uniref90_diamond "$HUMANN_DB_DIR"
humann_config --update database_folders nucleotide "$HUMANN_DB_DIR/chocophlan"
humann_config --update database_folders protein "$HUMANN_DB_DIR/uniref"
echo "HUMAnN Databases: $HUMANN_DB_DIR"

echo "Download eggNOG-mapper Data" | tee -a "$LOG_FILE"
mkdir -p "$EGGNOG_DIR"
"$CONDA_PREFIX/bin/download_eggnog_data.py" -y --dbname Bacteria -d 2 -P -M --data_dir "$EGGNOG_DIR"
export EGGNOG_DATA_DIR="$EGGNOG_DIR"
echo 'export EGGNOG_DATA_DIR="'$EGGNOG_DIR'"' >> ~/.bashrc

if [ -f "$EGGNOG_DATA_DIR/bacteria.dmnd" ]; then
  echo "Bacteria-only DIAMOND DB created at: $EGGNOG_DATA_DIR/bacteria.dmnd" | tee -a "$LOG_FILE"
else
  echo "ERROR: Expected bacteria.dmnd not found in $EGGNOG_DATA_DIR" | tee -a "$LOG_FILE"
  exit 1
fi

echo "eggNOG-mapper Database: $EGGNOG_DIR" | tee -a "$LOG_FILE"


echo "Install InterProScan v${INTERPROSCAN_VERSION} Manually" | tee -a "$LOG_FILE"
cd "$INSTALL_DIR"
mkdir -p "interproscan_download"
cd "interproscan_download"

wget -q "${INTERPROSCAN_URL}/interproscan-${INTERPROSCAN_VERSION}-64-bit.tar.gz"
wget -q "${INTERPROSCAN_URL}/interproscan-${INTERPROSCAN_VERSION}-64-bit.tar.gz.md5"

echo "Verifying checksum..."
md5sum -c "interproscan-${INTERPROSCAN_VERSION}-64-bit.tar.gz.md5"

echo "Extracting InterProScan..."
tar -xzf "interproscan-${INTERPROSCAN_VERSION}-64-bit.tar.gz"
mv "interproscan-${INTERPROSCAN_VERSION}" "$INTERPROSCAN_DIR"
export PATH="$INTERPROSCAN_DIR:$PATH"
echo 'export PATH="'$INTERPROSCAN_DIR':$PATH"' >> ~/.bashrc

echo "InterProScan version:" | tee -a "$LOG_FILE"
"$INTERPROSCAN_DIR/interproscan.sh" -version | tee -a "$LOG_FILE"


if [ -f "$INTERPROSCAN_DIR/setup.py" ]; then
  echo "Running InterProScan setup.py with interproscan.properties..." | tee -a "$LOG_FILE"
  python3 "$INTERPROSCAN_DIR/setup.py" -f "$INTERPROSCAN_DIR/interproscan.properties"
else
  echo "setup.py not found in $INTERPROSCAN_DIR — skipping this step." | tee -a "$LOG_FILE"
fi
echo "InterProScan: $INTERPROSCAN_DIR"


echo "Version Checks" | tee -a "$LOG_FILE"
{
  echo "FastQC: $(fastqc --version | head -n1)"
  echo "Cutadapt: $(cutadapt --version)"
  echo "TrimGalore: $(trim_galore --version)"
  echo "SortMeRNA: $(sortmerna --version)"
  echo "Bowtie2: $(bowtie2 --version)"
  echo "DIAMOND: $(diamond version | head -n1)"
  echo "eggNOG-mapper: $($CONDA_PREFIX/bin/emapper.py --version)"
  echo "MMseqs2: $(mmseqs --version | head -n1)"
  echo "MultiQC: $(multiqc --version)"
  echo "HMMER: $(hmmsearch -h 2>&1 | grep HMMER | head -n1)"
  echo "Prodigal: $(prodigal -v 2>&1 | grep -oE 'Prodigal V[0-9.]+' || echo 'not found')"
  echo "HUMAnN: $(humann --version)"
  echo "OpenJDK: $(java -version 2>&1 | head -n1)"
  echo "Python: $(python --version)"
  echo "Pandas: $(python -c 'import pandas; print(pandas.__version__)')"
  echo "NumPy: $(python -c 'import numpy; print(numpy.__version__)')"
  echo "SciPy: $(python -c 'import scipy; print(scipy.__version__)')"
  echo "Statsmodels: $(python -c 'import statsmodels; print(statsmodels.__version__)')"
  echo "scikit-learn: $(python -c 'import sklearn; print(sklearn.__version__)')"
  echo "Matplotlib: $(python -c 'import matplotlib; print(matplotlib.__version__)')"
  echo "Seaborn: $(python -c 'import seaborn; print(seaborn.__version__)')"
  echo "Biopython: $(python -c 'import Bio; print(Bio.__version__)')"
  echo "tqdm: $(python -c 'import tqdm; print(tqdm.__version__)')"
  echo "Requests: $(python -c 'import requests; print(requests.__version__)')"
  echo "lxml: $(python -c 'import lxml; print(lxml.__version__)')"
  echo "OpenPyXL: $(python -c 'import openpyxl; print(openpyxl.__version__)')"
} | tee -a "$LOG_FILE"

echo "Environment: $ENV_NAME"
echo "SETUP COMPLETE" | tee -a "$LOG_FILE"
