import argparse
import os
import time
from typing import Any, Dict, List

import pandas as pd
from Bio import Entrez

# NCBI requires email to work
Entrez.email = "omukutirodney@gmail.com"


def download_virus_sequences(virus: str, metadata_path: str, output_dir: str) -> None:
    """Downloads sequences for a given virus.

    Parameters
    ----------
    virus :  str
        Name of the virus.
    output_dir : str
        Path to where the sequences will be saved.
    """

    # Search NCBI accession ids for the sequences
    virus = virus.strip()
    search_handle = Entrez.esearch(db="nuccore", term=f"{virus}[orgn]", retmax=1000_000,
                                   idtype="acc")
    search_records = Entrez.read(search_handle)

    # Filter ids for whole genomes
    sequence_ids = list(filter(lambda x: x.startswith("NC_"), search_records["IdList"]))
    print(f"Downloading {virus} | Found {len(sequence_ids)} sequences")

    if len(sequence_ids):
        # Download sequences
        download_handle = Entrez.efetch(db="nucleotide", id=sequence_ids, retmode="xml")
        download_records = Entrez.read(download_handle)
        countries = list(map(extract_country, download_records))

        # Save downloaded sequences from the records
        for sequence_id, record in zip(sequence_ids, download_records):
            save_sequence(sequence_id, record, output_dir)

        # Save metadata
        save_metadata(sequence_ids, countries, download_records, metadata_path)


def extract_country(record: Dict[str, Any]) -> str:
    # Handle cases where the country does not exist
    try:
        features = record['GBSeq_feature-table'][0]['GBFeature_quals']
        features_dict = {}
        for feature in features:
            features_dict[feature['GBQualifier_name']] = feature['GBQualifier_value']

        return features_dict['country']
    except:
        return ""


def save_sequence(sequence_id: str, record: Dict[str, Any], output_dir: str) -> None:
    # Create directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    with open(os.path.join(output_dir, f"{sequence_id}.fa"), "w") as f:
        # Write sequence header
        f.write(f">{sequence_id} {record['GBSeq_definition']}\n")
        # Write the actual sequence
        f.write(record['GBSeq_sequence'].upper())


def save_metadata(sequence_ids: List[str], countries: List[str],
                  records: List[Dict[str, Any]], metadata_path: str) -> None:
    df = pd.DataFrame()
    # Add metadata to dataframe
    for sequence_id, country, record in zip(sequence_ids, countries, records):
        metadata_df = pd.DataFrame([{
            "id": sequence_id,
            "country": country,
            "moltype": record["GBSeq_moltype"],
            "source": record["GBSeq_source"],
            "taxonomy": record["GBSeq_taxonomy"]
        }])
        df = pd.concat([df, metadata_df], ignore_index=True)

    # Append to existing csv file
    df.to_csv(metadata_path, mode="a", index=False, header=False)


def read_file(filename: str) -> List[str]:
    with open(filename, 'r') as f:
        return f.readlines()

# download the sequences with metadata
def download_sequences(viruses_filepath: str, metadata_path: str,
                       data_dir: str) -> None:
    for virus in read_file(viruses_filepath):
        download_virus_sequences(virus, metadata_path, data_dir)
        time.sleep(5)

# specify the paths for your results
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("viruses_path", type=str,
                        help="path to txt file with virus names")
    parser.add_argument("metadata_path", type=str,
                        help="path where metadata will be saved")
    parser.add_argument("data_dir", type=str,
                        help="folder where data/sequences will be saved")
    args = parser.parse_args()

    download_sequences(args.viruses_path, args.metadata_path, args.data_dir)
