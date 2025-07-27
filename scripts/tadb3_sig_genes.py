import pandas as pd
from Bio import Entrez
from tqdm import tqdm
import re
from collections import Counter

Entrez.email = "Youremail@email.com" # Replace this
sig_path = "DA_results/sig_genes.txt"
tadb_path = "tadb3/master_table_tadb3.txt"
output_path = "tadb3/sig_genes_annotations.csv"

# Load significant UniRef90 IDs
sig_genes = set(pd.read_csv(sig_path, header=None)[0].str.strip())

# Load and filter TADB3 data
df = pd.read_csv(tadb_path, sep="\t")
df.columns = df.columns.str.strip()
df = df[df["UniRef90_ID"].isin(sig_genes)].copy()

# Ensure Score is numeric and filter by identity ≥ 90%
df["Bit_Score"] = pd.to_numeric(df["Bit_Score"], errors="coerce")
df = df[df["Bit_Score"] >= 90]
print(f"Total rows in master table: {len(df)}")
print(f"Significant UniRef90_IDs: {len(sig_genes)}")
print(f"Rows after UniRef90_ID filter: {df.shape[0]}")

# Extract species-level info from Organism field
df["Organism_str"] = df["Organism"]
df["Organism"] = df["Organism"].apply(lambda org: " ".join(str(org).split()[:2]))


# Fetch and clean protein descriptions from NCBI Entrez
def clean_protein_description(desc):
    desc = str(desc)
    desc = re.sub(r"\b(MULTISPECIES:|uncharacterized protein|hypothetical protein|putative)\b", "", desc, flags=re.IGNORECASE)
    desc = re.sub(r"\[.*?\]", "", desc)  # remove bracketed strain info
    return desc.strip().strip('.')

def fetch_protein_descriptions(accessions):
    acc_info = {}
    for acc in tqdm(accessions, desc="Fetching NCBI protein descriptions"):
        try:
            handle = Entrez.efetch(db="protein", id=acc, rettype="gb", retmode="text")
            record = handle.read()
            desc = ""
            for line in record.splitlines():
                if line.startswith("DEFINITION"):
                    desc = line.replace("DEFINITION", "").strip().rstrip(".")
                    break
            acc_info[acc] = clean_protein_description(desc)
        except Exception:
            acc_info[acc] = ""
    return acc_info

df["Protein_Accession"] = df["Protein_Accession"].astype(str)
desc_dict = fetch_protein_descriptions(df["Protein_Accession"].unique())
df["Protein_Description"] = df["Protein_Accession"].map(lambda x: desc_dict.get(x, ""))

# Aggregation helper functions
def custom_agg(x, mode="collapse"):
    if mode == "collapse":
        return "; ".join(sorted(set(map(str, x.dropna()))))
    elif mode == "mean":
        return x.astype(float).mean()
    elif mode == "min":
        return x.astype(float).min()
    elif mode == "max":
        return x.astype(float).max()
    elif mode == "first":
        return str(x.dropna().iloc[0]) if not x.dropna().empty else ""

def most_common_entry(x):
    values = [str(i).strip() for i in x.dropna()]
    return Counter(values).most_common(1)[0][0] if values else ""

def prioritize_validation(x):
    return "Experimental" if "Experimental" in set(x) else "Computational"

# Aggregation map
agg_funcs = {
    "Query_ID": lambda x: custom_agg(x),
    "Bit_Score": lambda x: custom_agg(x, "max"),
    "Alignment_Length": lambda x: custom_agg(x, "mean"),
    "Mismatch": lambda x: custom_agg(x, "mean"),
    "Gap_Open": lambda x: custom_agg(x, "first"),
    "Q_Start": lambda x: custom_agg(x, "first"),
    "Q_End": lambda x: custom_agg(x, "first"),
    "S_Start": lambda x: custom_agg(x, "first"),
    "S_End": lambda x: custom_agg(x, "first"),
    "E_Value": lambda x: custom_agg(x, "min"),
    "Sample_Name": lambda x: custom_agg(x),
    "Validation_Type": prioritize_validation,
    "Toxin_Antitoxin": lambda x: custom_agg(x),
    "Dataset": lambda x: custom_agg(x),
    "Protein_Accession": lambda x: custom_agg(x),
    "Genome_Accession": lambda x: custom_agg(x),
    "Genome_Position": lambda x: custom_agg(x),
    "Organism": lambda x: custom_agg(x),
    "Organism_str": lambda x: custom_agg(x),
    "Protein_Description": most_common_entry,
}

# Collapse by UniRef90_ID
collapsed_df = df.groupby("UniRef90_ID", as_index=False).agg(agg_funcs)

# Count distinct samples
collapsed_df["Sample_Count"] = collapsed_df["Sample_Name"].apply(lambda x: len(set(str(x).split(";"))))

# Save to CSV
collapsed_df.to_csv(output_path, index=False)
print(f"\nTADB3 annotations table saved to:\n{output_path}")
