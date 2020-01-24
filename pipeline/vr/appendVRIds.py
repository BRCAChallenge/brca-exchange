#!/usr/bin/env python
import requests
import logging
import sys
import argparse
import csv
# from ga4gh.core import sha512t24u, ga4gh_digest, ga4gh_identify, ga4gh_serialize
# from ga4gh.vr import __version__, models, normalize
# from ga4gh.vr.extras.dataproxy import SeqRepoRESTDataProxy
# from ga4gh.vr.extras.translator import Translator


# SEQREPO_REST_SERVICE_URL = "http://localhost:5000/seqrepo"
# DP = SeqRepoRESTDataProxy(base_url=SEQREPO_REST_SERVICE_URL)
# TLR = Translator(data_proxy=DP,
#                  translate_sequence_identifiers=True,
#                  normalize=True,
#                  identify=True)


def parse_args():
    parser = argparse.ArgumentParser(description='Determine correct VR id.')
    parser.add_argument('-l', '--logfile', default='/tmp/append_vr_ids.log')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'),
                        help='Input variants.')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                        help='Output variants.')
    options = parser.parse_args()
    return options


def main(args):
    options = parse_args()
    inputFile = options.input
    outputFile = options.output
    logfile = options.logfile

    logging.basicConfig(filename=logfile, filemode="w", level=logging.WARNING,
                        format=' %(asctime)s %(filename)-15s %(message)s')
    input_file = csv.reader(inputFile, delimiter='\t')
    output_file = csv.writer(outputFile, delimiter='\t')
    
    for row in input_file:
        input_header_row = row
        break

    new_column_to_append = ["VR_ID"]

    output_header_row = input_header_row + new_column_to_append

    output_file.writerow(output_header_row)

    hgvsIndex = input_header_row.index("Genomic_HGVS_38")
    
    for variant in input_file:
        hgvs = variant[hgvsIndex]
        
        # Add empty data by default
        variant.append('-')

        if not is_empty(hgvs):
            vr_id = get_vr_id(hgvs)
            variant[output_header_row.index("VR_ID")] = vr_id

        output_file.writerow(variant)


def is_empty(field_value):
    return field_value == '' or field_value is None or field_value == '-'

def get_vr_id(hgvs):
    try:
        # allele = TLR.from_hgvs(hgvs)
        # allele_dict = allele.as_dict()
        # return allele_dict['_id']
        return "placeholder"
    except ValueError as e:
        logging.warning("Exception during processing of " + str(hgvs) + ": " + str(e))
        return '-'


if __name__ == "__main__":
    sys.exit(main(sys.argv))
