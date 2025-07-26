import os
import glob
import pandas as pd
import re
from tqdm import tqdm
from Bio import Entrez

 
Entrez.email = "shrivishalinir@gmail.com"
sig_genes_path = "tadb3/sig_genes.txt"
vfdb_paths = [
    "vfdb/dieguez/*_vs_VFDB.tsv",
    "vfdb/ev/*_vs_VFDB.tsv"
]
output_path = "vfdb/sig_genes_annotations.csv"

# Load significant UniRef90_IDs 
sig_genes = set(pd.read_csv(sig_genes_path, header=None)[0].str.strip())

# Fetch gene/protein info from NCBI 
def fetch_ncbi_info(gb_ids):
    gb_info = {}
    for gb_id in tqdm(gb_ids, desc="Fetching NCBI gb| info"):
        try:
            handle = Entrez.efetch(db="protein", id=gb_id, rettype="gb", retmode="text")
            record = handle.read()
            gene = ""
            desc = ""
            for line in record.splitlines():
                if line.startswith("DEFINITION"):
                    desc = line.replace("DEFINITION", "").strip().rstrip(".")
                if "/gene=" in line:
                    gene = line.split('=')[-1].strip('"')
            gb_info[gb_id] = {"Gene_Name": gene, "Protein_Description": desc}
        except Exception:
            gb_info[gb_id] = {"Gene_Name": "", "Protein_Description": ""}
    return gb_info

# Load and filter all VFDB records 
vfdb_records = []

for path in vfdb_paths:
    for file in glob.glob(path):
        dataset = "dieguez" if "dieguez" in file else "ev"
        sample = os.path.basename(file).split("_vs_VFDB.tsv")[0]

        df = pd.read_csv(file, sep="\t", header=None)
        df.columns = ['UniRef90_ID', 'Query_ID', 'Percent_Identity', 'Alignment_Length', 'Mismatch', 'Gap_Open',
                      'Q_Start', 'Q_End', 'S_Start', 'S_End', 'E_Value', 'Bit_Score']

        df = df[df["UniRef90_ID"].isin(sig_genes)].copy()
        if df.empty:
            continue

        # Extract gb accession
        df["gb_accession"] = df["Query_ID"].apply(lambda x: re.search(r"gb\|([\w_\.]+)", x).group(1) if "gb|" in x else "")
        df["Sample_Name"] = sample
        df["Dataset"] = dataset

        vfdb_records.append(df)

if not vfdb_records:
    print("No matching records found.")
    exit()

vfdb_df = pd.concat(vfdb_records, ignore_index=True)

# Fetch gene/protein annotations 
vfdb_df["gb_accession"] = vfdb_df["gb_accession"].astype(str)
unique_gb = vfdb_df["gb_accession"].unique()
gb_info = fetch_ncbi_info(unique_gb)

vfdb_df["Gene_Name"] = vfdb_df["gb_accession"].map(lambda x: gb_info.get(x, {}).get("Gene_Name", ""))
vfdb_df["Protein_Description"] = vfdb_df["gb_accession"].map(lambda x: gb_info.get(x, {}).get("Protein_Description", ""))

# Collapse and aggregate 
def custom_agg(x, mode="collapse"):
    if mode == "collapse":
        return "; ".join(sorted(set(map(str, x.dropna()))))
    elif mode == "mean":
        return x.astype(float).mean()
    elif mode == "min":
        return x.astype(float).min()
    elif mode == "max":
        return x.astype(float).max()

agg_funcs = {
    "Query_ID": lambda x: custom_agg(x, "collapse"),
    "Percent_Identity": lambda x: custom_agg(x, "max"),
    "Alignment_Length": lambda x: custom_agg(x, "mean"),
    "Mismatch": lambda x: custom_agg(x, "mean"),
    "Gap_Open": lambda x: custom_agg(x, "collapse"),
    "Q_Start": lambda x: custom_agg(x, "collapse"),
    "Q_End": lambda x: custom_agg(x, "collapse"),
    "S_Start": lambda x: custom_agg(x, "collapse"),
    "S_End": lambda x: custom_agg(x, "collapse"),
    "E_Value": lambda x: custom_agg(x, "min"),
    "Bit_Score": lambda x: custom_agg(x, "max"),
    "Sample_Name": lambda x: custom_agg(x, "collapse"),
    "Dataset": lambda x: custom_agg(x, "collapse"),
    "Gene_Name": lambda x: custom_agg(x, "collapse"),
    "Protein_Description": lambda x: custom_agg(x, "collapse"),
}

vfdb_collapsed = vfdb_df.groupby("UniRef90_ID", as_index=False).agg(agg_funcs)
vfdb_collapsed["Sample_Count"] = vfdb_collapsed["Sample_Name"].apply(lambda x: len(set(str(x).split(";"))))

# Save final table 
vfdb_collapsed.to_csv(output_path, index=False)
