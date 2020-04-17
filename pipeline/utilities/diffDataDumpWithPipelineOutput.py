#!/usr/bin/env python

'''
This script compares CSV output from the DB to a tsv file (e.g. built_with_change_types.tsv or built.tsv).

To use, follow these steps:

1. Create a CSV file from the DB (simply update the "Data_Release_id" to correspond to a specific data release, leave
   the fields the same):

    `Copy (SELECT "id", "Genomic_Coordinate_hg38", "Pathogenicity_all", "Pathogenicity_expert", "Source",
    "HGVS_cDNA" from variant where "Data_Release_id"=1 ORDER BY "id" ASC) To '/tmp/pgdump.csv' CSV HEADER;`

2. Update the headers in the pgdump.csv (dump from DB) file to match their respective headers in the built.tsv
   (pipeline output) file. e.g. Genomic_Coordinate_hg38 maps to pyhgvs_Genomic_Coordinate_38 in the tsv file output
   from the pipeline.

3. Double check that all necessary fields are present or adjust the script as needed!!

4. Run `python diffDataDumpWithPipelineOutput.py --built /PATH/TO/built.tsv --pgdump /PATH/TO/pgdump.csv`
'''

import argparse
import csv


def addGsIfNecessary(genomic_coordinate):
    if ":g." not in genomic_coordinate:
            genomic_coordinate = genomic_coordinate[:6] + 'g.' + genomic_coordinate[6:]
    return genomic_coordinate


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--built", default="built.tsv",
                        help="Built file from pipeline output")
    parser.add_argument("--pgdump", default="pgdump.csv",
                        help="Dump of PG data in CSV format")

    args = parser.parse_args()
    built = csv.DictReader(open(args.built, "r"), delimiter="\t")
    pgDump = csv.DictReader(open(args.pgdump, "r"), delimiter=",")

    bad_matches = {}

    built_variants = {}
    pg_variants = {}

    for built_variant in built:
        built_genomic_coordinate = addGsIfNecessary(built_variant["pyhgvs_Genomic_Coordinate_38"])
        built_pathogenicity_all = built_variant["Pathogenicity_all"]
        built_pathogenicity_expert = built_variant["Pathogenicity_expert"]
        built_cDNA = built_variant["pyhgvs_cDNA"]
        built_sources = built_variant["Source"]
        built_variants[built_genomic_coordinate] = (built_pathogenicity_all, built_pathogenicity_expert, built_cDNA, built_sources)

    for pg_variant in pgDump:
        pg_genomic_coordinate = pg_variant["pyhgvs_Genomic_Coordinate_38"]
        pg_pathogenicity_all = pg_variant["Pathogenicity_all"]
        pg_pathogenicity_expert = pg_variant["Pathogenicity_expert"]
        pg_cDNA = pg_variant["pyhgvs_cDNA"]
        pg_sources = pg_variant["Source"]
        pg_variants[pg_genomic_coordinate] = (pg_pathogenicity_all, pg_pathogenicity_expert, pg_cDNA, pg_sources)

    for variant in list(built_variants.keys()):
        if built_variants[variant] != pg_variants[variant]:
            print("BAD MATCH for variant %s: %s | %s" % (variant, built_variants[variant], pg_variants[variant]))
            bad_matches[variant] = (built_variants[variant], pg_variants[variant])

    if bad_matches == {}:
        print("All variants in pgDump match!")


if __name__ == "__main__":
    main()
