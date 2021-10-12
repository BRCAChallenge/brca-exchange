#!/usr/bin/env python
import requests
import logging
import sys
import argparse
import csv
import os
import vcf
from copy import deepcopy
from ga4gh.core import sha512t24u, ga4gh_digest, ga4gh_identify, ga4gh_serialize
from ga4gh.vrs import __version__, models, normalize
from ga4gh.vrs.dataproxy import SeqRepoRESTDataProxy
from ga4gh.vrs.extras.translator import Translator

csv.field_size_limit(10000000)

SEQREPO_REST_SERVICE_URL = "http://localhost:5000/seqrepo"
DP = SeqRepoRESTDataProxy(base_url=SEQREPO_REST_SERVICE_URL)
TLR = Translator(data_proxy=DP,
                 translate_sequence_identifiers=True,
                 normalize=True,
                 identify=True)

def parse_args():
    parser = argparse.ArgumentParser(description='Determine correct VR id.')
    parser.add_argument('-l', '--logfile', default='/tmp/append_vr_ids.log')
    parser.add_argument('-i', '--inputdir', help='Input file directory',
                        default="/home/brca/pipeline-data/pipeline-input/")
    parser.add_argument('-a', '--artifacts', help='Artifacts directory',
                        default="/home/brca/pipeline-data/artifacts/"),
    parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                        help='Output variants.')
    options = parser.parse_args()
    return options


def main(args):
    options = parse_args()
    input_dir = options.inputdir
    output_file = options.output
    artifacts_dir = options.artifacts
    logfile = options.logfile

    logging.basicConfig(filename=logfile, filemode="w", level=logging.WARNING,
                        format=' %(asctime)s %(filename)-15s %(message)s')

    output_file = csv.writer(output_file, delimiter='\t')
    print(f'input dir={input_dir}')
    print(f'outputfile={output_file}')
    print(f'art dir={artifacts_dir}')
    input_files = []
    # get files
    for filename in os.listdir(input_dir):
        # TODO: add support for tsv files
        if filename.endswith(".vcf"):
            input_files.append(os.path.join(input_dir, filename))
        else:
            continue

    submissions = {}

    for file in input_files:
        # identify source
        filename = file.split('/')[-1]
        source = get_source(filename)
        print(source)
        # ensure each entry is an individual submission
        f_in = open(file, "r")
        f_out = open(os.path.join(artifacts_dir, source + "ready.vcf"), "w")
        one_variant_transform(f_in, f_out, source)
        f_in.close()
        f_out.close()
        print('done w transform')
        # get individual submissions to retrieve vrs id
        f = open(os.path.join(artifacts_dir, source + "ready.vcf"), "r")
        vcf_reader = vcf.Reader(f, strict_whitespace=True)
        for record in vcf_reader:
            print(record)
            # get info necessary to get VRS id
            ref = record.REF.replace("-", "")
            v = f'{record.CHROM}-{record.POS}-{ref}-{record.ALT[0]}'
            # get vrs id
            if not is_empty(v):
                print(v)
                vrs_id = get_vrs_id(v)
                print(vrs_id)
                # variant[output_header_row.index("VR_ID")] = vrs_id
            # format fields
            # write submission to output_file
        sys.exit(1)

def one_variant_transform(f_in, f_out, source_name):
    """takes a vcf file, read each row, if the ALT field contains more than
       one item, create multiple variant row based on that row. also adds
       ids to all individual reports (each line in the vcf). writes new vcf"""
    vcf_reader = vcf.Reader(f_in, strict_whitespace=True)
    vcf_writer = vcf.Writer(f_out, vcf_reader)
    count = 1
    for record in vcf_reader:
        n = len(record.ALT)
        if n == 1:
            if source_name == "ExAC":
                record = append_exac_allele_frequencies(record)
            record.INFO['BX_ID'] = count
            count += 1
            vcf_writer.write_record(record)
        else:
            for i in range(n):
                new_record = deepcopy(record)
                new_record.ALT = [deepcopy(record.ALT[i])]
                new_record.INFO['BX_ID'] = count
                count += 1
                for key in record.INFO.keys():
                    value = deepcopy(record.INFO[key])
                    if type(value) == list and len(value) == n:
                        new_record.INFO[key] = [value[i]]
                if source_name == "ExAC":
                    new_record = append_exac_allele_frequencies(record, new_record, i)
                vcf_writer.write_record(new_record)


def append_exac_allele_frequencies(record, new_record=None, i=None):
    if new_record is None:
        for subpopulation in EXAC_SUBPOPULATIONS:
            # calculate allele frequencies for each subpopulation
            allele_count = record.INFO[("AC_" + subpopulation)]
            allele_number = record.INFO[("AN_" + subpopulation)]
            allele_frequency = "-"
            if len(allele_count) > 0 and allele_number != 0:
                allele_frequency = float(allele_count[0]) / float(allele_number)
                allele_frequency = str(utilities.round_sigfigs(allele_frequency, 3))
            record.INFO[("AF_" + subpopulation)] = allele_frequency
        return record
    else:
        new_record.INFO['AF'] = record.INFO['AF'][i]
        for subpopulation in EXAC_SUBPOPULATIONS:
            allele_count = record.INFO[("AC_" + subpopulation)][i]
            allele_number = record.INFO[("AN_" + subpopulation)]
            allele_frequency = "-"
            if allele_number != 0:
                allele_frequency = float(allele_count) / float(allele_number)
                allele_frequency = str(utilities.round_sigfigs(allele_frequency, 3))
            new_record.INFO[("AF_" + subpopulation)] = allele_frequency
        return new_record


def get_source(filename):
    # TODO: implement
    return filename.split('.')[0]

def is_empty(field_value):
    return field_value == '' or field_value is None or field_value == '-'

def get_vrs_id(v):
    # TODO: test
    a = TLR.translate_from("NC_000013.11:g.32315668G>C","hgvs")
    print(a.as_dict())
    sys.exit(1)
    print(translate_from("NC_000013.11:g.32936732G>C","hgvs"))
    try:
        allele = TLR.translate_from(v, "gnomad")
        print(f"allele: {allele}")
        allele_dict = allele.as_dict()
        print(f"alleledict: {allele_dict}")
        return allele_dict['_id']
    except IndexError as e:
        print(e)
        logging.warning("Exception during processing of " + str(v) + ": " + str(e))
        return '-'


if __name__ == "__main__":
    sys.exit(main(sys.argv))
