"""
Adding CA_ID from clingen AlleleRegistry API (http://reg.clinicalgenome.org/doc/AlleleRegistry_1.01.xx_api_v1.pdf)
"""
import argparse
import logging
import re
import sys
import time
from typing import List, Iterable

import pandas as pd
import requests

from common import utils

HGVS_COL = 'Genomic_HGVS_38'
CA_ID_COL = 'CA_ID'

ALLELES_ENDPOINT = "http://reg.clinicalgenome.org/alleles?file=hgvs"
MAX_TRIES = 5
CA_ID_RE = re.compile(r'CA[0-9]+')

RESP_ID_KEY = '@id'
RESP_ERROR_TYPE_KEY = 'errorType'
RESP_DESCRIPTION_KEY = 'description'
RESP_MESSAGE_KEY = 'message'

def parse_args():
    parser = argparse.ArgumentParser(description='Gathers ClinGen Allele Registry IDs.')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'),
                        help='Input variants.')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                        help='Output variants.')
    parser.add_argument('-l', '--log', type=str,
                        help='Log file')
    options = parser.parse_args()
    return options


def _extract_ca_id(r, orig_hgvs):
    if RESP_ID_KEY not in r.keys():
        if RESP_ERROR_TYPE_KEY in r.keys():
            logging.warning(
                f"Server could not process {orig_hgvs}: {r[RESP_ERROR_TYPE_KEY]}. {r[RESP_DESCRIPTION_KEY]}. {r[RESP_MESSAGE_KEY]}")
        return None

    id_str = r[RESP_ID_KEY]

    ca_id = id_str.split('/')[-1]

    if CA_ID_RE.match(ca_id):
        return ca_id
    else:
        logging.warning(f"Could not extract ca_id for {orig_hgvs} out of {ca_id}")
        return None


def _perform_request(url, data):
    tries = 0

    req = None
    while not req and tries < MAX_TRIES:
        try:
            req = requests.post(url, data)
        except requests.exceptions.RequestException as e:
            logging.error(f"Request to {url} failed. Going to retry")
            logging.exception(e)
            time.sleep(10)
            tries += 1

    if not req:
        logging.error(f"Requests failed {MAX_TRIES} times, exiting.")
        sys.exit(1)

    if req.status_code != 200:
        logging.error(f"Request failed, server responded with status code {req.status_code} and message {req.json()}")
        sys.exit(1)

    return req.json()


def _query_server(hgvs_ids: Iterable[str], url: str, max_chunk_size: int = 50000) -> List[str]:
    resp_all = []
    for split in utils.split_list_in_chunks(hgvs_ids, max_chunk_size):
        hgvs_data = '\n'.join(split)
        logging.info(f"Requesting chunk of size {len(split)} from {url}")
        resp_all.extend(_perform_request(url, hgvs_data))

    return [_extract_ca_id(r, h) for r, h in zip(resp_all, hgvs_ids)]


def add_caid_to_build_df(df: pd.DataFrame, url: str) -> pd.DataFrame:
    hgvs_avail = ~df[HGVS_COL].isna()

    df.loc[hgvs_avail, CA_ID_COL] = _query_server(df.loc[hgvs_avail, HGVS_COL].values, url)

    return df


def main(args):
    options = parse_args()
    input_file = options.input
    output_file = options.output

    utils.setup_logfile(options.log)

    input_df = pd.read_csv(input_file, sep='\t')
    output_df = add_caid_to_build_df(input_df, ALLELES_ENDPOINT)

    output_df.to_csv(output_file, sep='\t', index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
