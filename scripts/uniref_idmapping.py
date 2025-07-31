import os
import argparse
import requests
from datetime import datetime


def extract_uniref90_ids(input_dir, output_dir):
    ids_dir = os.path.join(output_dir, "extracted_ids")
    os.makedirs(ids_dir, exist_ok=True)

    for filename in os.listdir(input_dir):
        if filename.endswith("_genefamilies.tsv"):
            sample_name = filename.replace("_genefamilies.tsv", "")
            ids_file = os.path.join(ids_dir, f"{sample_name}_uniref90_ids.txt")
            unique_ids = set()

            with open(os.path.join(input_dir, filename), "r") as file:
                for line in file:
                    if line.startswith("UniRef90"):
                        unique_ids.add(line.split("\t")[0].split("|")[0])  # Extract UniRef90 ID

            with open(ids_file, "w") as f:
                f.write("\n".join(sorted(unique_ids)))

            print(f"Extracted {len(unique_ids)} UniRef90 IDs for {sample_name} to {ids_file}")


def split_ids_into_chunks(ids, chunk_size):
    for i in range(0, len(ids), chunk_size):
        yield ids[i:i + chunk_size]


def download_fasta_for_sample(ids_file, fasta_dir, log_dir, chunk_size=500):
    sample_name = os.path.basename(ids_file).replace("_uniref90_ids.txt", "")
    fasta_file = os.path.join(fasta_dir, f"{sample_name}.fasta")
    log_file = os.path.join(log_dir, f"{sample_name}.log")

    with open(ids_file, "r") as f:
        uniref90_ids = f.read().splitlines()

    with open(log_file, "w") as log:
        log.write(f"Start processing {sample_name} at {datetime.now()}\n")
        log.write(f"Total UniRef90 IDs: {len(uniref90_ids)}\n")

        for i, chunk in enumerate(split_ids_into_chunks(uniref90_ids, chunk_size)):
            query = " OR ".join(chunk)
            base_url = "https://rest.uniprot.org/uniref/stream"
            params = {"format": "fasta", "query": query}

            log.write(f"Processing chunk {i + 1} with {len(chunk)} IDs...\n")
            response = requests.get(base_url, params=params)

            if response.status_code == 200:
                with open(fasta_file, "a") as fasta_out:
                    fasta_out.write(response.text)
                log.write(f"Chunk {i + 1} completed successfully.\n")
            else:
                log.write(f"Chunk {i + 1} failed. Status Code: {response.status_code}\n")

        log.write(f"Finished processing {sample_name} at {datetime.now()}\n")
    print(f"FASTA sequences saved to {fasta_file}, log written to {log_file}")


def main():
    parser = argparse.ArgumentParser(description="Extract and download UniRef90 FASTA sequences per file.")
    parser.add_argument("--input_dir", required=True, help="Directory containing Humann3 genefamilies TSV files.")
    parser.add_argument("--output_dir", required=True, help="Directory to save extracted IDs, FASTA files, and logs.")
    parser.add_argument("--chunk_size", type=int, default=500, help="Number of IDs per API query.")
    args = parser.parse_args()

    # Create output subdirectories
    fasta_dir = os.path.join(args.output_dir, "fastas")
    log_dir = os.path.join(args.output_dir, "logs")
    os.makedirs(fasta_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    # Extract UniRef90 IDs for each file
    extract_uniref90_ids(args.input_dir, args.output_dir)

    # Download FASTA sequences for each file
    ids_dir = os.path.join(args.output_dir, "extracted_ids")
    for ids_file in os.listdir(ids_dir):
        if ids_file.endswith("_uniref90_ids.txt"):
            download_fasta_for_sample(os.path.join(ids_dir, ids_file), fasta_dir, log_dir, args.chunk_size)


if __name__ == "__main__":
    main()
