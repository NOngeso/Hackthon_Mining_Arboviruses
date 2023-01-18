import argparse
from typing import Any, Dict, List

import pandas as pd
from Bio import Entrez, SeqIO

# retry failed requests
Entrez.max_tries = 10


class MetadataDownloader:
    def __init__(self, viruses_path: str, countries_path: str, metadata_path: str,
                 email: str = None, api_key: str = None,
                 retmax: int = 1_000_000) -> None:
        self.viruses_path = viruses_path
        self.countries_path = countries_path
        self.metadata_path = metadata_path
        self.retmax = retmax

        # setup NCBI access requirements
        Entrez.email = email
        Entrez.api_key = api_key

        self.metadata_df = pd.DataFrame()

    def run(self) -> None:
        try:
            # read a list countries and viruses
            countries = MetadataDownloader.read_lines(self.countries_path)
            viruses = MetadataDownloader.read_lines(self.viruses_path)

            # download metadata for each virus at a time for a given country
            for country in countries:
                for virus in viruses:
                    self.__fetch_by_virus_and_country(virus, country)

        except:
            print('Failed to download metadata')

        # save downloaded metadata
        self.metadata_df.set_index('id', inplace=True)
        self.metadata_df.to_csv(self.metadata_path)

    @staticmethod
    def read_lines(filename: str) -> List[str]:
        with open(filename, 'r') as f:
            # reads lines without the newline character
            return f.read().splitlines()

    def __fetch_by_virus_and_country(self, virus: str, country: str) -> None:
        # fetch sequences ids for the virus
        sequences_ids = self.__fetch_sequences_ids(virus, country)
        # fetch records for the viruses sequences_ids
        if len(sequences_ids):
            print(
                f'Fetching metadata for: {virus} from: {country} | Found {len(sequences_ids)} sequences'
            )
            try:
                records = self.__fetch_records(sequences_ids)
                self.__extract_metadata(records)
            except:
                print(f'Failed to fetch and extract metadata for: {virus}')
        else:
            print(f'No results for: {virus} from: {country}')

    def __fetch_sequences_ids(self, virus: str, country: str) -> List[str]:
        handle = Entrez.esearch(db='nucleotide',
                                term=f'{virus}[orgn] AND {country}[country]',
                                retmax=self.retmax)
        return Entrez.read(handle)['IdList']

    def __fetch_records(self, sequences_ids: List[str]) -> List[Dict[str, Any]]:
        records = []
        for id in sequences_ids:
            handle = Entrez.efetch(db='nucleotide', id=id, rettype='gb')
            records.append(SeqIO.read(handle, 'gb'))

        return records

    def __extract_metadata(self, records: List[SeqIO.SeqRecord]) -> None:
        for record in records:
            features = {k: v[0] for k, v in record.features[0].qualifiers.items()}
            df = pd.DataFrame([{
                'id': record.id,
                'description': record.description,
                'pubmed_id': record.annotations['references'][0].pubmed_id,
                'sequence_length': len(record.seq),
                **features
            }])

            self.metadata_df = pd.concat([self.metadata_df, df], ignore_index=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('viruses_path', type=str,
                        help='path to txt file with virus names separated by a line')
    parser.add_argument('countries_path', type=str,
                        help='path to txt file with country names separated by a line')
    parser.add_argument('metadata_path', type=str,
                        help='path to where metadata will be saved')
    parser.add_argument('--email', type=str, help='email to give to NCBI', default=None)
    parser.add_argument('--api_key', type=str, help='api key for NCBI', default=None)
    parser.add_argument('--retmax', type=int, help='maximum number of results',
                        default=1_000_000)
    args = parser.parse_args()

    metadata_downloader = MetadataDownloader(args.viruses_path, args.countries_path,
                                             args.metadata_path, args.email,
                                             args.api_key, args.retmax)
    metadata_downloader.run()
