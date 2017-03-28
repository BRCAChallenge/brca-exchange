#!/usr/bin/env python
"""
checks that all reports are accounted for in final pipeline output
"""
import argparse
import csv
import logging
from os import listdir
from os.path import isfile, join, abspath
from aggregate_reports import get_reports_files
import vcf


ARGS = None


def main():
    global ARGS
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--built",
                        help="built.tsv with all variants")
    parser.add_argument("-r", "--ready_input_dir",
                        help="file directory with all procesed files with bx_ids used to compile built.tsv")
    parser.add_argument('-a', "--artifacts_dir", help='Artifacts directory with pipeline artifact files.')
    parser.add_argument("-v", "--verbose", action="count", default=True, help="determines logging")

    ARGS = parser.parse_args()

    configure_logging()

    bx_ids = get_bx_ids()

    matches_per_source = find_matches_per_source(bx_ids)

    missing_reports = find_missing_reports(matches_per_source, bx_ids)

    logging.debug("Reports absent from release: %s", missing_reports)


def get_bx_ids():
    # Get all bx_ids present in source files organized by source
    bx_ids = {}

    files = get_reports_files(ARGS.ready_input_dir)

    for file in files:
        file_path = abspath(ARGS.ready_input_dir + file)
        if file_path.endswith('.tsv'):
            source = "ENIGMA"
            bx_ids[source] = []
            tsv_file = csv.DictReader(open(file_path, "r"), delimiter='\t')
            for report in tsv_file:
                ids = map(int, report['BX_ID'].split(','))
                bx_ids[source] = bx_ids[source] + ids
        else:
            suffix = '.vcf'
            source = file[:(len(file)-len(suffix))]
            bx_ids[source] = []
            vcf_reader = vcf.Reader(open(file_path, 'r'), strict_whitespace=True)
            try:
                for record in vcf_reader:
                    ids = map(int, record.INFO['BX_ID'])
                    bx_ids[source] = bx_ids[source] + ids
            except ValueError as e:
                print e

    return bx_ids


def find_matches_per_source(bx_ids):
    matches_per_source = {}
    built = csv.DictReader(open(ARGS.built, "r"), delimiter='\t')
    column_prefix = "BX_ID_"
    bx_id_columns = [f for f in built.fieldnames if column_prefix in f]
    matches_per_source = {}
    for variant in built:
        for column in bx_id_columns:
            source = column[len(column_prefix):]
            if source not in matches_per_source:
                matches_per_source[source] = []
            source_bx_ids = variant[column]
            variant_sources = variant["Source"].split(',')
            match = False
            for src in variant_sources:
                if source == src:
                    match = True
            if isEmpty(source_bx_ids):
                if match:
                    logging.warning("Variant %s has source %s but no report ids from that source", variant, source)
            else:
                if not match:
                    logging.warning("Variant %s has report(s) %s from source %s, but source is not associated with variant", variant, source_bx_ids, source)
                else:
                    source_bx_ids = map(int, source_bx_ids.split(','))
                    for source_bx_id in source_bx_ids:
                        if source_bx_id in bx_ids[source]:
                            matches_per_source[source].append(source_bx_id)
                        else:
                            logging.warning("Report(s) %s found on variant %s, but report does not exist from source %s", source_bx_ids, variant, source)
    return matches_per_source


def find_missing_reports(matches_per_source, bx_ids):
    missing_reports = {}
    for source in matches_per_source:
        matches_set = get_set_of_ids(matches_per_source[source])
        original_ids_set = get_set_of_ids(bx_ids[source])
        missing_reports[source] = matches_set.symmetric_difference(original_ids_set)
    return missing_reports


def get_set_of_ids(ids):
    id_set = set(ids)
    assert len(ids) == len(set(id_set))
    return id_set


def configure_logging():
    if ARGS.verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.CRITICAL

    log_file_path = ARGS.artifacts_dir + "missing_reports.log"
    logging.basicConfig(filename=log_file_path, filemode="w", level=logging_level)


def isEmpty(value):
    return value == '-' or value is None or value == '' or value == []


if __name__ == "__main__":
    main()
