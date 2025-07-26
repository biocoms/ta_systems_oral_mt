import os
import glob
import pandas as pd
from collections import defaultdict, Counter
from tqdm import tqdm
from io import StringIO

 
sig_genes_path = "DA_results/sig_genes.txt"
eggnog_dirs = {
    "dieguez": "eggnog/dieguez/*_eggnog.emapper.annotations",
    "ev": "eggnog/ev/*_eggnog.emapper.annotations"
}
output_csv = "eggnog/sig_genes_annotations.csv"

# Load significant UniRef90 IDs 
sig_genes = set(
    pd.read_csv(sig_genes_path, header=None)[0].str.strip().str.split('.').str[0]
)

# Function to correctly parse EggNOG-mapper file 
def parse_emapper(file):
    with open(file) as f:
        lines = [line for line in f if not line.startswith("##")]

    # Find the true header line (starts with "#query")
    header_index = next(i for i, line in enumerate(lines) if line.startswith("#query"))
    header = lines[header_index].strip().lstrip('#').split('\t')
    data_lines = lines[header_index + 1:]

    from io import StringIO
    df = pd.read_csv(StringIO(''.join(data_lines)), sep='\t', names=header)
    return df


# Dictionary to collect annotations 
annotated_genes = defaultdict(list)

# Parse EggNOG files 
for dataset, pattern in eggnog_dirs.items():
    for file in tqdm(glob.glob(pattern), desc=f"Reading {dataset} files"):
        sample = os.path.basename(file).split('_eggnog')[0]
        try:
            df = parse_emapper(file)

            if "query" not in df.columns:
                print(f"[ERROR] 'query' column missing in file: {file}")
                print(f"Columns found: {df.columns.tolist()}")
                continue

            df = df[df["query"].str.split('.').str[0].isin(sig_genes)].copy()
            df["Dataset"] = dataset
            df["Sample"] = sample

            for _, row in df.iterrows():
                uid_clean = str(row["query"]).split(".")[0]
                annotated_genes[uid_clean].append(row)

        except Exception as e:
            print(f"[ERROR] Failed to read {file}: {e}")
            continue

# Collapse annotations per UniRef90 ID 
collapsed_data = []

for uid, records in annotated_genes.items():
    samples = {r["Sample"] for r in records}
    datasets = {r["Dataset"] for r in records}
    sample_count = len(samples)

    df_records = pd.DataFrame(records)
    row_out = {
        "UniRef90ID": uid,
        "Dataset": ";".join(sorted(datasets)),
        "Samples": ";".join(sorted(samples)),
        "Sample_Count": sample_count
    }

    for col in df_records.columns:
        if col in ["query", "Dataset", "Sample"]:
            continue
        values = df_records[col]
        values = values[values != "-"]
        if not values.empty:
            row_out[col] = Counter(values).most_common(1)[0][0]
        else:
            row_out[col] = "-"

    collapsed_data.append(row_out)

# Save final DataFrame 
final_df = pd.DataFrame(collapsed_data)
final_df.to_csv(output_csv, index=False)
print(f"\n EggNOG annotations table saved to: {output_csv}")
