import os
import argparse

import pandas as pd
from Bio import Entrez

# retry failed requests
Entrez.max_tries = 10


class GBDownloader:
    def __init__(self, metadata_path: str, output_dir: str, email: str = None,
                 api_key: str = None, retmax: int = 1_000_000) -> None:
        self.metadata_path = metadata_path
        self.output_dir = output_dir
        self.retmax = retmax

        # setup NCBI access requirements
        Entrez.email = email
        Entrez.api_key = api_key

        # Load sequence metadata
        metadata_df = pd.read_csv(metadata_path)

        # Filter complete genomes based on description
        complete_genomes_df = metadata_df[metadata_df['description'].str.contains(
            'complete genome')]

        self.sequences_ids = complete_genomes_df['id']

    def run(self) -> None:
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)

        for id in self.sequences_ids:
            handle = Entrez.efetch(db='nucleotide', id=id, rettype='gb')
            with open(os.path.join(self.output_dir, f'{id}.gb'), 'w') as f:
                f.write(handle.read())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('metadata_path', type=str,
                        help='path to where metadata will be saved')
    parser.add_argument('output_dir', type=str, help='path to where data will be saved')
    parser.add_argument('--email', type=str, help='email to give to NCBI', default=None)
    parser.add_argument('--api_key', type=str, help='api key for NCBI', default=None)
    parser.add_argument('--retmax', type=int, help='maximum number of results',
                        default=1_000_000)
    args = parser.parse_args()

    gb_downloader = GBDownloader(args.metadata_path, args.output_dir, args.email,
                                 args.api_key, args.retmax)
    gb_downloader.run()
