import os
import pandas as pd
from Bio import SeqIO, Entrez
from tqdm import tqdm


Entrez.email = "shrivishalinir@gmail.com"  # Replace this
diamond_dirs = {"dieguez": "tadb3/diamond_dieguez", "ev": "tadb3/diamond_ev"}
fasta_file = "tadb3/tadb3_combined.fasta"
sig_ids_file = "DA_results/sig_genes.txt"
output_master_prefix = "master_table"
output_filtered = "tadb3/tadb3_tox_antitox_sig_with_protein_names.tsv"
output_main = "tadb3/tadb3_tox_antitox_main_table.tsv"

# Parse FASTA metadata 
def parse_fasta(fasta_file):
    records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        parts = record.description.split("|")
        if len(parts) >= 5:
            records.append({
                "Query_ID": record.id,
                "Protein_Accession": parts[0],
                "Genome_Accession": parts[1],
                "Genome_Position": parts[2],
                "Organism": parts[4]
            })
    return pd.DataFrame(records)

# Build per-dataset master tables 
def build_master_tables(diamond_dirs, fasta_df):
    dfs = []
    for dataset, path in diamond_dirs.items():
        for file in os.listdir(path):
            if not file.endswith(".tsv"):
                continue
            fpath = os.path.join(path, file)
            sample = file.replace(".tsv", "")
            df = pd.read_csv(fpath, sep="\t", header=None)
            df.columns = ["Query_ID", "UniRef90_ID", "Identity", "Alignment_Length", "Mismatch",
                          "Gap_Open", "Q_Start", "Q_End", "S_Start", "S_End", "E_Value", "Bit_Score"]
            df = df.merge(fasta_df, on="Query_ID", how="left")
            df["Sample_Name"] = sample
            df["Dataset"] = dataset
            df["Validation_Type"] = df["Query_ID"].apply(lambda x: "Experimental" if "exp" in x else "Computational")
            df["Toxin_Antitoxin"] = df["Query_ID"].apply(lambda x: "Toxin" if "_tox_" in x else "Antitoxin")
            dfs.append(df)
    return pd.concat(dfs)

# Write master tables 
def save_master_tables(df, prefix):
    for dataset in df["Dataset"].unique():
        output_path = f"tadb3/{prefix}_{dataset}.tsv"
        df[df["Dataset"] == dataset].to_csv(output_path, sep="\t", index=False)
    
    # Save combined master table
    df.to_csv("tadb3/master_table_tadb3.txt", sep="\t", index=False)


# Filter significant genes 
with open(sig_ids_file, "r") as f:
    sig_ids = set(line.strip() for line in f if line.strip())

# Fetch protein names using Entrez 
def fetch_ncbi_name(accession):
    try:
        handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
        for line in handle:
            if line.startswith("DEFINITION"):
                return line.replace("DEFINITION  ", "").strip().rstrip(".")
    except:
        return None

def annotate_protein_names(df):
    accs = df["Protein_Accession"].dropna().unique()
    name_map = {}
    for acc in tqdm(accs, desc="Fetching protein names"):
        name_map[acc] = fetch_ncbi_name(acc)
    df["Protein_Name"] = df["Protein_Accession"].map(name_map)
    return df

# Summarize into main table 
def summarize(df):
    summary = []
    for uid in df["UniRef90_ID"].unique():
        sub = df[df["UniRef90_ID"] == uid]
        name = sub["Protein_Name"].mode().iloc[0] if not sub["Protein_Name"].isnull().all() else "Unknown"
        species = "; ".join(sorted(sub["Organism"].dropna().unique()))
        bit_score = sub["Bit_Score"].max()
        identity = sub["Identity"].max()
        validation = "Both" if set(sub["Validation_Type"]) == {"Experimental", "Computational"} else sub["Validation_Type"].iloc[0]
        datasets = "; ".join(sorted(sub["Dataset"].unique()))
        nsamples = sub["Sample_Name"].nunique()
        summary.append({
            "UniRef90_ID": uid,
            "Protein_Name": name,
            "Representative_Species": species,
            "Best_Bit_Score": bit_score,
            "Best_Identity": identity,
            "Validation_Type": validation,
            "Dataset": datasets,
            "Unique_Samples": nsamples
        })
    return pd.DataFrame(summary)


def run_pipeline():
    fasta_df = parse_fasta(fasta_file)
    print("Parsed FASTA.")
    
    master_df = build_master_tables(diamond_dirs, fasta_df)
    print("Built master tables.")
    
    save_master_tables(master_df, output_master_prefix)
    print("Saved master tables.")
    
    filtered = master_df[master_df["UniRef90_ID"].isin(sig_ids)].copy()
    print(f"Filtered {len(filtered)} hits.")
    
    filtered_named = annotate_protein_names(filtered)
    filtered_named.to_csv(output_filtered, sep="\t", index=False)
    print("Saved filtered table with protein names.")
    
    final = summarize(filtered_named)
    final.to_csv(output_main, sep="\t", index=False)
    print("Main summary table written.")


run_pipeline()
