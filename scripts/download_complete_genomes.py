import argparse
import os
from typing import Any, Dict

import pandas as pd
from Bio import Entrez


def extract_sequence(sequence_id: str, record: Dict[str, Any]) -> str:
    record_str = f">{sequence_id} {record['GBSeq_definition']}\n"
    record_str += f"{record['GBSeq_sequence'].upper()}\n"

    return record_str


def save(data: str, filepath: str) -> None:
    with open(os.path.join(filepath), 'w') as f:
        f.write(data)


def download_complete_genomes(metadata_path: str, output_dir: str) -> None:
    # Load sequence metadata
    metadata_df = pd.read_csv(metadata_path)

    # Filter complete genomes based on description
    complete_genomes_df = metadata_df[metadata_df['description'].str.contains(
        'complete genome')]
    print(f'Found {len(complete_genomes_df)} complete genomes')
    print(complete_genomes_df.groupby(['organism'])['organism'].count())

    # Fetch all sequences at once
    sequences_ids = complete_genomes_df['id']
    download_handle = Entrez.efetch(db="nucleotide", id=sequences_ids, retmode="xml")
    sequences_records = Entrez.read(download_handle)

    # Save downloaded sequences as a combined file
    combined_sequences = ''.join(
        map(lambda x: extract_sequence(x[0], x[1]), zip(sequences_ids,
                                                        sequences_records)))
    # Create directory to hold the combined sequences
    combined_output_dir = os.path.join(output_dir, 'combined')
    os.makedirs(combined_output_dir, exist_ok=True)

    save(combined_sequences, os.path.join(combined_output_dir, 'combined_sequences.fa'))

    # Create directory to hold the separate sequences
    sequences_output_dir = os.path.join(output_dir, 'fasta')
    os.makedirs(sequences_output_dir, exist_ok=True)
    # Save downloaded sequences in separate files
    for sequence_id, record in zip(sequences_ids, sequences_records):
        output_path = os.path.join(sequences_output_dir, f'{sequence_id}.fa')
        save(extract_sequence(sequence_id, record), output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('metadata_path', type=str,
                        help='path to the csv file with metadata')
    parser.add_argument(
        'output_dir', type=str,
        help='path to directory that will contain the downloaded sequences')
    parser.add_argument('--email', type=str, help='entrez requires an email',
                        default='')

    args = parser.parse_args()
    Entrez.email = args.email
    download_complete_genomes(args.metadata_path, args.output_dir)