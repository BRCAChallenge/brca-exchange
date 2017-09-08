#!/usr/bin/env python
"""
remove field from tsv file
"""
import argparse
import csv

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        help="input file with merged ENIGMA data")
    parser.add_argument("-o", "--output",
                        help="Output file with corrected ENIGMA data")
    parser.add_argument("-f", "--field",
                        help="Field to remove")
    args = parser.parse_args()

    with open(args.input, 'r') as fin, open(args.output, 'w') as fout:
        reader = csv.reader(fin, dialect='excel-tab')
        writer = csv.writer(fout, dialect='excel-tab')
        index_of_field_to_remove = None
        for row in reader:
            if not index_of_field_to_remove:
                index_of_field_to_remove = row.index(args.field)
            del row[index_of_field_to_remove:index_of_field_to_remove + 1]
            writer.writerow(row)


if __name__ == "__main__":
    main()
