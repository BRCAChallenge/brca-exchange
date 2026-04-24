#!/usr/bin/env python
import argparse
import logging
import pysam
import os
import re

from common import config
from data_merging import utilities

from data_merging.variant_merging_constants import (
    DEFAULT_CONTENTS,
    FIELD_DICT,
    ENIGMA_FILE,
    COLUMN_SOURCE,
    VARIANT_REPORT_FILES
)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", help="Pathname of the config file")
    parser.add_argument("-i", "--input_dir", help="Directory with input VCF files")
    parser.add_argument("-l", "--log_file", help="Pathname of the log file", default="aggregate_reports.log")
    parser.add_argument("-m", "--merged_file", help="Pathname of the file with merged variants")
    parser.add_argument("-o", "--output_file", help="Output filename")
    parser.add_argument("-v", "--verbose", action="count", default=True, help="determines logging")
    args = parser.parse_args()
    return(args)
    

def main():
    args = parse_args()
    if args.verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.CRITICAL
    logging.basicConfig(filename=args.log_file, filemode="w", level=logging_level,
                        format=' %(asctime)s %(filename)-15s %(message)s')
    gene_config_df = config.load_config(args.config)
    genome_regions_symbol_dict = config.get_genome_regions_symbol_dict(gene_config_df, 'start_hg38_legacy_variants', 'end_hg38_legacy_variants')
    #
    # Read the first line to get the list of columns
    with open(args.merged_file, "r") as fp:
        line = fp.readline().strip()
        columns = re.split('\t', line)
    #
    # Process subsequent lines
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
        logging.info("finished normalizing %s with %d reports" % (file, len(reports)))
    reports_output_fp.close()
    logging.info("Done")


def aggregate_reports(reports_files, columns, genome_regions_symbol_dict):
    all_reports = []
    for file in reports_files:
        reports = normalize_reports(file, columns, genome_regions_symbol_dict)
        all_reports.extend(reports)
    return all_reports


def get_reports_files(input_directory):
    reports_files = []
    for f in os.listdir(input_directory):
        filename, file_extension = os.path.splitext(f)
        if (filename in FIELD_DICT and file_extension == ".vcf") or f == ENIGMA_FILE:
            reports_files.append(f)
    return reports_files




def normalize_reports(file, columns, genome_regions_symbol_dict):
    if os.path.basename(file) == ENIGMA_FILE:
        reports = normalize_enigma_vcf_reports(file, columns)
    elif os.path.splitext(file)[1] == ".vcf":
        filename = os.path.splitext(os.path.basename(file))[0]
        reports = normalize_vcf_reports(file, columns, filename, ".vcf",
                                        genome_regions_symbol_dict)
    else:
        raise Exception("ERROR: unexpected file format: %s" % file)
    return reports


def normalize_vcf_reports(file, columns, filename, file_extension, genome_regions_symbol_dict):
    reports = []
    reader = pysam.VariantFile(file)
    count = 0
    source_suffix = ".vcf"
    basename = os.path.basename(file)[:-len(source_suffix)]
    source = basename.split(".")[0]
    for record in reader:
        count += 1
        genome_coor = ("chr" + str(record.chrom) + ":g." + str(record.pos) + ":" +
                       record.ref + ">" + str(record.alts[0]))

        if utilities.is_outside_boundaries(record.chrom, record.pos, genome_regions_symbol_dict):
            logging.warning("Skipping report since the positions is outside the genome boundaries: " + str(record))
            continue
        report = utilities.associate_chr_pos_ref_alt_with_item(record, len(columns), source, genome_coor, genome_regions_symbol_dict)
        for key, value in FIELD_DICT[source].items():
            try:
                column_name = key + "_" + source
                column_index = columns.index(column_name)
                report[column_index] = record.info[value]
            except KeyError:
                if value in reader.header.info:
                    # Field declared in header but absent from this record (e.g. gnomAD
                    # FAF fields not computed for every variant).
                    column_name = key + "_" + source
                    column_index = columns.index(column_name)
                    report[column_index] = DEFAULT_CONTENTS
                else:
                    raise Exception("WARNING: Key error with report: %s \n\nError on value: %s \n\n Error in record.info: %s \n\nNeeds attn." % (report, value, record.info))
        reports.append(report)
    reader.close()
    return reports


def normalize_enigma_vcf_reports(file, columns):
    reports = []
    vcf_reader = pysam.VariantFile(file)
    info_field_names = list(vcf_reader.header.info.keys())
    enigma_columns = utilities.add_columns_to_enigma_data_vcf(info_field_names)
    enigma_column_indexes = {key: value for key, value in enumerate(enigma_columns)}
    for record in vcf_reader:
        (items, chrom, pos, ref, alt) = \
            utilities.associate_chr_pos_ref_alt_with_enigma_vcf_record(
                record, info_field_names)
        report = [DEFAULT_CONTENTS] * len(columns)
        for key, value in enigma_column_indexes.items():
            if value in columns:
                report[columns.index(value)] = items[key]
        reports.append(report)
    vcf_reader.close()
    return reports


def write_report(report, columns, reports_output_fp):
    if len(report) != len(columns):
        raise Exception("mismatching number of columns in head and row")
    for ii in range(len(report)):
        val = report[ii]
        if isinstance(val, (list, tuple)):
            report[ii] = ",".join(str(xx) for xx in val)
        elif not isinstance(val, str):
            report[ii] = DEFAULT_CONTENTS if val is None else str(val)
    reports_output_fp.write("\t".join(report)+"\n")
    

if __name__ == "__main__":
    main()
