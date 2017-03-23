#!/usr/bin/env python
"""
this scripts takes the enigma variant list and merge vcf files in a folder into
the exisitng enigma variants:
"""
import argparse
import os
import vcf
import logging
import csv
from variant_merging import GENOME1K_FIELDS, CLINVAR_FIELDS, LOVD_FIELDS, EX_LOVD_FIELDS, BIC_FIELDS, ESP_FIELDS, EXAC_FIELDS, FIELD_DICT, ENIGMA_FILE, COLUMN_SOURCE, COLUMN_GENE, COLUMN_GENOMIC_HGVS, COLUMN_VCF_CHR, COLUMN_VCF_POS, COLUMN_VCF_REF, COLUMN_VCF_ALT


def write_reports_tsv(filename, columns, ready_files_dir):
    reports_output = open(filename, "w")
    reports_output.write("\t".join(columns)+"\n")

    reports_files = [ready_files_dir + r for r in get_reports_files(ready_files_dir)]

    reports = aggregate_reports(reports_files, columns)

    for report in reports:
        if len(reports) != len(columns):
            raise Exception("mismatching number of columns in head and row")
        for ii in range(len(reports)):
            if type(reports[ii]) == list:
                comma_delimited_string = ",".join(str(xx) for xx in reports[ii])
                reports[ii] = comma_delimited_string
            elif type(reports[ii]) == int:
                reports[ii] = str(reports[ii])
        reports_output.write("\t".join(reports)+"\n")
    reports_output.close()

    print "final number of reports: %d" % len(reports)
    print "Done"


def aggregate_reports(reports_files, columns):
    # Gathers all reports from an input directory, normalizes them, and combines them into a single list.

    for file in reports_files:
        file_reports = normalize_reports(file, columns)
        print "finished normalizing %s" % (file)
        reports = reports + file_reports

    return reports


def get_reports_files(input_directory):
    return [f for f in os.listdir(input_directory) if os.path.isfile(os.path.join(input_directory, f)) and "ready" in f or "ENIGMA_combined_with_bx_ids" in f]


def normalize_reports(file, columns):
    reports = []
    filename, file_extension = os.path.splitext(file)
    # with open(file, "r") as f_in:
    if file_extension == ".vcf":
        reader = vcf.Reader(open(file, "r"), strict_whitespace=True)
        count = 0
        source_suffix = "ready.vcf"
        source = os.path.basename(file)[:-len(source_suffix)]
        for record in reader:
            count += 1
            genome_coor = ("chr" + str(record.CHROM) + ":g." + str(record.POS) + ":" +
                           record.REF + ">" + str(record.ALT[0]))

            # first set all values to default
            report = ['-'] * len(columns)

            report[COLUMN_SOURCE] = source
            if record.CHROM == "13":
                report[COLUMN_GENE] = "BRCA2"
            elif record.CHROM == "17":
                report[COLUMN_GENE] = "BRCA1"
            else:
                raise Exception("Wrong chromosome")
            report[COLUMN_GENOMIC_HGVS] = genome_coor
            report[COLUMN_VCF_CHR] = record.CHROM
            report[COLUMN_VCF_POS] = record.POS
            report[COLUMN_VCF_REF] = record.REF
            report[COLUMN_VCF_ALT] = str(record.ALT[0])
            for value in FIELD_DICT[source].values():
                try:
                    report.append(record.INFO[value])
                except KeyError:
                    # logging.warning("KeyError appending VCF record.INFO[value] to variant. Variant: %s \n Record.INFO: %s \n value: %s", variants[genome_coor], record.INFO, value)
                    # if source == "BIC":
                        # variants[genome_coor].append(DEFAULT_CONTENTS)
                        # logging.debug("Could not find value %s for source %s in variant %s, inserting default content %s instead.", value, source, DEFAULT_CONTENTS)
                    # else:
                        # raise Exception("There was a problem appending a value for %s to variant %s" % (value, variants[genome_coor]))
                    print "WARNING: Key error, needs attn."
            reports.append(report)
    # TODO:
    # elif file_extension == ".tsv":
        # handle enigma data
    return reports


if __name__ == "__main__":
    main()
