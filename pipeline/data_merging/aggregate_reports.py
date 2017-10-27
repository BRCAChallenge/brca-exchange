#!/usr/bin/env python
import argparse
import os
import vcf
import logging
import csv
import urllib
from variant_merging import (
     add_columns_to_enigma_data,
     associate_chr_pos_ref_alt_with_enigma_item,
     associate_chr_pos_ref_alt_with_item,
     BIC_FIELDS,
     COLUMN_SOURCE,
     COLUMN_GENE,
     COLUMN_GENOMIC_HGVS,
     COLUMN_VCF_CHR,
     COLUMN_VCF_POS,
     COLUMN_VCF_REF,
     COLUMN_VCF_ALT,
     CLINVAR_FIELDS,
     EX_LOVD_FIELDS,
     ESP_FIELDS,
     EXAC_FIELDS,
     FIELD_DICT,
     ENIGMA_FILE,
     GENOME1K_FIELDS,
     LOVD_FIELDS
     )


def write_reports_tsv(filename, columns, ready_files_dir):
    reports_output = open(filename, "w")

    reports_files = [ready_files_dir + r for r in get_reports_files(ready_files_dir)]

    reports = aggregate_reports(reports_files, columns)

    reports_output.write("\t".join(columns)+"\n")

    for report in reports:
        if len(report) != len(columns):
            raise Exception("mismatching number of columns in head and row")
        for ii in range(len(report)):
            if type(report[ii]) == list:
                comma_delimited_string = ",".join(str(xx) for xx in report[ii])
                report[ii] = comma_delimited_string
            elif type(report[ii]) == int:
                report[ii] = str(report[ii])
        reports_output.write("\t".join(report)+"\n")

    reports_output.close()

    print "final number of reports: %d" % len(reports)
    print "Done"


def aggregate_reports(reports_files, columns):
    # Gathers all reports from an input directory, normalizes them, and combines them into a single list.
    reports = []

    for file in reports_files:
        file_reports = normalize_reports(file, columns)
        print "finished normalizing %s" % (file)
        reports = reports + file_reports

    return reports


def get_reports_files(input_directory):
    reports_files = []
    for f in os.listdir(input_directory):
        filename, file_extension = os.path.splitext(f)
        if (filename in FIELD_DICT and file_extension == ".vcf") or f == ENIGMA_FILE:
            reports_files.append(f)
    return reports_files


def normalize_reports(file, columns):
    filename, file_extension = os.path.splitext(file)
    if file_extension == ".vcf":
        reports = normalize_vcf_reports(file, columns, filename, file_extension)
    elif file_extension == ".tsv":
        if os.path.basename(file) != ENIGMA_FILE:
            raise Exception("ERROR: received tsv file that is not for ENIGMA: %s" % (file))
        reports = normalize_enigma_tsv_reports(file, columns, filename, file_extension)
    for report in reports:
        if len(report) != len(columns):
            report += [DEFAULT_CONTENTS] * len(FIELD_DICT[source])
    for report in reports:
        if len(report) != len(columns):
            raise Exception("mismatching number of columns in head and row")
    return reports


def normalize_vcf_reports(file, columns, filename, file_extension):
    reports = []
    if "clinvar" in filename.lower():
        # Descriptions and Summary Evidence in clinvar contain spaces -- cause strict whitespace failure
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

        report = associate_chr_pos_ref_alt_with_item(record, len(columns), source, genome_coor)
        for key, value in FIELD_DICT[source].iteritems():
            try:
                column_name = key + "_" + source
                column_index = columns.index(column_name)
                if source == "LOVD":
                    report[column_index] = map(urllib.unquote_plus, record.INFO[value])
                else:
                    report[column_index] = record.INFO[value]
            except KeyError:
                raise Exception("WARNING: Key error with report: %s \n\nError on value: %s \n\n Error in record.INFO: %s \n\nNeeds attn." % (report, value, record.INFO))
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
            enigma_columns = add_columns_to_enigma_data(line)
            for key, value in enumerate(enigma_columns):
                enigma_column_indexes[key] = value
        else:
            (items, chrom, pos, ref, alt) = associate_chr_pos_ref_alt_with_enigma_item(line)
            report = ['-'] * len(columns)
            for key, value in enigma_column_indexes.iteritems():
                report[columns.index(value)] = items[key]
            reports.append(report)
    enigma_file.close()
    return reports
