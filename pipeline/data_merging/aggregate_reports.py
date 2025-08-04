#!/usr/bin/env python
import argparse
import logging
import os
import re
import vcf

from common import config
from data_merging import utilities

from data_merging.variant_merging_constants import (
    DEFAULT_CONTENTS,
    FIELD_DICT,
    ENIGMA_FILE,
    COLUMN_SOURCE
)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", help="Directory with input VCF files")
    parser.add_argument("-o", "--output_file", help="Output filename")
    parser.add_argument("-c", "--config", help="Pathname of the config file")
    parser.add_argument("-m", "--merged_file", help="Pathname of the file with merged variants")
    parser.add_argument("-v", "--verbose", action="count", default=False, help="determines logging")
    args = parser.parse_args()
    return(args)
    

def main():
    args = parse_args()
    gene_config_df = config.load_config(args.config)
    genome_regions_symbol_dict = config.get_genome_regions_symbol_dict(gene_config_df, 'start_hg38_legacy_variants', 'end_hg38_legacy_variants')
    with open(args.merged_file, "r") as fp:
        line = fp.readline().strip()
        columns = re.split('\t', line)
    write_reports_tsv(args.output_file, columns, args.input_dir, genome_regions_symbol_dict)
    

def write_reports_tsv(filename, columns, ready_files_dir,
                      genome_regions_symbol_dict):
    # For each data source, normalize the contents of the reports, the data submissions,
    # so that for each variant, each report is represented by a consistent set of columns.
    reports_output_fp = open(filename, "w")
    reports_output_fp.write("\t".join(columns)+"\n")
    reports_files = [ready_files_dir + r for r in get_reports_files(ready_files_dir)]
    for file in reports_files:
        source_file = os.path.basename(file)
        source_name = re.split("\.", source_file)
        logging.info("Starting to normalize reports for %s (%s)" % (source_name, file))
        reports = normalize_reports(file, columns, genome_regions_symbol_dict)
        for this_report in reports:
            # If any columns are missing, pad with the default value.
            # This isn't needed with ENIGMA data, which has a fixed TSV format.
            if source_file != ENIGMA_FILE:
                if len(this_report) != len(columns):
                    this_report += [DEFAULT_CONTENTS] * len(FIELD_DICT[source])
            write_report(this_report, columns, reports_output_fp)
        logging.info("finished normalizing %s" % (file))
    reports_output_fp.close()
    logging.info("final number of reports: %d" % len(reports))
    logging.info("Done")


def get_reports_files(input_directory):
    reports_files = []
    for f in os.listdir(input_directory):
        filename, file_extension = os.path.splitext(f)
        if (filename in FIELD_DICT and file_extension == ".vcf") or f == ENIGMA_FILE:
            reports_files.append(f)
    return reports_files


def normalize_reports(file, columns, genome_regions_symbol_dict):
    filename, file_extension = os.path.splitext(file)
    if file_extension == ".vcf":
        reports = normalize_vcf_reports(file, columns, filename, file_extension, genome_regions_symbol_dict)
    elif file_extension == ".tsv":
        if os.path.basename(file) != ENIGMA_FILE:
            raise Exception("ERROR: received tsv file that is not for ENIGMA: %s" % (file))
        reports = normalize_enigma_tsv_reports(file, columns, filename, file_extension)
    return reports


def normalize_vcf_reports(file, columns, filename, file_extension, genome_regions_symbol_dict):
    reports = []
    if "clinvar" in filename.lower():
        # If fields contain spaces they cause strict whitespace failure
        strict_whitespace = False
    else:
        strict_whitespace = True
    reader = vcf.Reader(open(file, "r"), strict_whitespace=strict_whitespace)
    count = 0
    source_suffix = ".vcf"
    source = os.path.basename(file)[:-len(source_suffix)]
    for record in reader:
        count += 1
        genome_coor = ("chr" + str(record.CHROM) + ":g." + str(record.POS) + ":" +
                       record.REF + ">" + str(record.ALT[0]))
        if utilities.is_outside_boundaries(record.CHROM, record.POS, genome_regions_symbol_dict):
            logging.warning("Skipping report since the positions is outside the genome boundaries: " + str(record))
            continue
        report = utilities.associate_chr_pos_ref_alt_with_item(record, len(columns), source, genome_coor, genome_regions_symbol_dict)
        for key, value in FIELD_DICT[source].items():
            try:
                column_name = key + "_" + source
                column_index = columns.index(column_name)
                report[column_index] = record.INFO[value]
            except KeyError:
                logging.warning("Key error with report: %s \n\nError on key %s, value: %s variant: %s\n"
                                % (report, key, value, genome_coor))
                report[column_index] = DEFAULT_CONTENTS
        logging.info("source %s record number %d coord %s" % (source, count, genome_coor))
        reports.append(report)
    return reports


def normalize_enigma_tsv_reports(file, columns, filename, file_extension):
    reports = []
    enigma_file = open(file, 'r')
    line_num = 0
    enigma_column_indexes = {}
    for line in enigma_file:
        line_num += 1
        if line_num == 1:
            enigma_columns = utilities.add_columns_to_enigma_data(line)
            for key, value in enumerate(enigma_columns):
                enigma_column_indexes[key] = value
        else:
            (items, chrom, pos, ref, alt) = utilities.associate_chr_pos_ref_alt_with_enigma_item(line)
            report = ['-'] * len(columns)
            for key, value in enigma_column_indexes.items():
                report[columns.index(value)] = items[key]
            reports.append(report)
    enigma_file.close()
    return reports


def write_report(report, columns, reports_output_fp):
    if len(report) != len(columns):
        raise Exception("mismatching number of columns in head and row")
    for ii in range(len(report)):
        if type(report[ii]) == list:
            comma_delimited_string = ",".join(str(xx) for xx in report[ii])
            report[ii] = comma_delimited_string
        elif type(report[ii]) == int:
            report[ii] = str(report[ii])
    reports_output_fp.write("\t".join(report)+"\n")
    

if __name__ == "__main__":
    main()
