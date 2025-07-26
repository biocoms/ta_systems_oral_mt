import os
import glob
import pandas as pd
import re
import requests
from tqdm import tqdm
from collections import defaultdict, Counter


sig_path = "DA_results/sig_genes.txt"
interpro_paths = [
    "interproscan/dieguez/*.tsv",
    "interproscan/ev/*.tsv"
]
output_csv = "interproscan/sig_genes_annotations.csv"


# Load significant UniRef90 IDs
sig_genes = set(pd.read_csv(sig_path, header=None)[0].str.strip())

# Read and combine InterProScan .tsv files
interpro_records = []
for path in interpro_paths:
    for file in glob.glob(path):
        dataset = "dieguez" if "dieguez" in file else "ev"
        sample = os.path.basename(file).replace(".tsv", "")
        df = pd.read_csv(file, sep="\t", header=None, dtype=str)
        df.columns = [
            "Protein_Accession", "Sequence_MD5", "Sequence_Length", "Analysis",
            "Signature_Accession", "Signature_Description", "Start", "End", "Score",
            "Status", "Date", "InterPro_Accession", "InterPro_Description", "GO", "Pathways"
        ]
        df["Dataset"] = dataset
        df["Sample_Name"] = sample
        interpro_records.append(df)

full_df = pd.concat(interpro_records, ignore_index=True)
full_df = full_df[full_df["Protein_Accession"].isin(sig_genes)].copy()

# Extract GO IDs
def extract_go_ids(go_str):
    if pd.isna(go_str): return []
    return re.findall(r"GO:\d{7}", go_str)

full_df["GO_IDs"] = full_df["GO"].apply(extract_go_ids)
all_go_ids = sorted(set(go for ids in full_df["GO_IDs"] for go in ids))

# Query GO terms from QuickGO
def query_quickgo(go_ids):
    go_info = {}
    url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
    headers = {"Accept": "application/json"}
    for go_id in tqdm(go_ids, desc="Querying GO terms"):
        try:
            r = requests.get(url + go_id, headers=headers)
            if r.ok:
                result = r.json()["results"][0]
                go_info[go_id] = {
                    "name": result.get("name", ""),
                    "aspect": result.get("aspect", "")
                }
        except:
            go_info[go_id] = {"name": "", "aspect": ""}
    return go_info

go_annotations = query_quickgo(all_go_ids)

# Annotate GO categories
def map_go_annotations(go_ids):
    all_terms, bp, mf, cc = [], [], [], []
    for go in go_ids:
        entry = go_annotations.get(go, {})
        name = entry.get("name", "")
        aspect = entry.get("aspect", "")
        if name: all_terms.append(name)
        if aspect == "biological_process": bp.append(name)
        elif aspect == "molecular_function": mf.append(name)
        elif aspect == "cellular_component": cc.append(name)
    return (
        "; ".join(sorted(set(all_terms))),
        "; ".join(sorted(set(bp))),
        "; ".join(sorted(set(mf))),
        "; ".join(sorted(set(cc)))
    )

full_df[["GO_Terms", "GO_BP", "GO_MF", "GO_CC"]] = full_df["GO_IDs"].apply(lambda ids: pd.Series(map_go_annotations(ids)))

# Collapse all information
grouped = defaultdict(dict)
for _, row in full_df.iterrows():
    uid = row["Protein_Accession"]
    analysis = row["Analysis"]

    # Analysis-specific fields
    grouped[uid].setdefault(f"{analysis}_acc", []).append(row["Signature_Accession"])
    grouped[uid].setdefault(f"{analysis}_desc", []).append(row["Signature_Description"])
    try:
        grouped[uid].setdefault(f"{analysis}_evalue", []).append(float(row["Score"]))
    except:
        grouped[uid].setdefault(f"{analysis}_evalue", []).append(1e6)

    # GO terms
    for col in ["GO_Terms", "GO_BP", "GO_MF", "GO_CC"]:
        grouped[uid].setdefault(col, []).extend(row[col].split("; ") if pd.notna(row[col]) else [])

    # MetaCyc/Reactome
    pathways = str(row["Pathways"]).split("|") if pd.notna(row["Pathways"]) else []


    # Additional metadata
    grouped[uid].setdefault("Dataset", set()).add(row["Dataset"])
    grouped[uid].setdefault("Sample_Name", set()).add(row["Sample_Name"])
    grouped[uid].setdefault("Sequence_Lengths", []).append(int(row["Sequence_Length"]))

# Assemble final table
final_rows = []
for uid, data in grouped.items():
    row = {
        "UniRef90_ID": uid,
        "Dataset": ";".join(sorted(data["Dataset"])),
        "Sample_Count": len(data["Sample_Name"]),
        "Samples": ";".join(sorted(data["Sample_Name"])),
        "Avg_Sequence_Length": int(sum(data["Sequence_Lengths"]) / len(data["Sequence_Lengths"])),
        "GO_Terms": "; ".join(sorted(set(data["GO_Terms"]))),
        "GO_BP": "; ".join(sorted(set(data["GO_BP"]))),
        "GO_MF": "; ".join(sorted(set(data["GO_MF"]))),
        "GO_CC": "; ".join(sorted(set(data["GO_CC"])))
    }
    for key in data:
        if key.endswith("_acc"):
            row[key] = "; ".join(sorted(set(data[key])))
        elif key.endswith("_desc"):
            row[key] = Counter(data[key]).most_common(1)[0][0]
        elif key.endswith("_evalue"):
            row[key] = min(data[key])  # Best (lowest) e-value
    final_rows.append(row)

# Save to CSV
final_df = pd.DataFrame(final_rows)
final_df.to_csv(output_csv, index=False)
print(f"InterProScan annotation tables saved to:\n{output_csv}")
