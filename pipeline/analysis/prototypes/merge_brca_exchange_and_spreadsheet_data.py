#!/usr/bin/env python3
"""
Merge V4_Popfreq data from built_with_v4b.tsv into ENIGMA_Johanna.tsv

This script:
1. Reads built_with_v4b.tsv and creates a lookup dictionary
2. For each row, takes column 405 (pyhgvs_cDNA) and replaces the part before
   the first colon with column 383 (Gene_Symbol) to create a join key
3. Maps this key to column 523 (V4_Popfreq)
4. Reads ENIGMA_Johanna.tsv and adds the V4_Popfreq column by matching
   column 1 (Variant) with the join key
"""

import argparse
import csv
import sys

import openpyxl
import xlrd


def get_allele_count_value(faf95_popmax_population, row):
    """
    Create allele count and allele number column indices from the faf95_popmax_population value
    and return the corresponding values from the row.

    Args:
        faf95_popmax_population: The population code (e.g., "nfe", "afr")
        row: Dictionary containing the row data

    Returns:
        Tuple of (allele_count_value, allele_number_value), or ('', '') if not found.
        All values default to empty strings.
    """
    if not faf95_popmax_population or faf95_popmax_population == '-' or faf95_popmax_population == 'None':
        return '', ''

    population_upper = faf95_popmax_population.upper()
    allele_count_index = f"Allele_count_joint_{population_upper}_GnomADv4"
    allele_number_index = f"Allele_number_joint_{population_upper}_GnomADv4"

    allele_count = row.get(allele_count_index, '')
    allele_number = row.get(allele_number_index, '')

    # Convert None to empty string
    return allele_count if allele_count else '', allele_number if allele_number else ''


def parse_enigma_class1_comment(comment):
    """
    Parse Comment_on_clinical_significance_ENIGMA for Class 1 frequency data.

    For comments beginning with "Class 1 not pathogenic based on frequency >1% in an outbred sampleset.",
    extracts the frequency and population from the second sentence (e.g. "Frequency 0.03659 (African)").

    Returns:
        Tuple of (frequency, population) strings, or (None, None) if the comment does not match.
    """
    class1_prefix = "Class 1 not pathogenic based on frequency >1% in an outbred sampleset."
    if not comment or not comment.startswith(class1_prefix):
        return None, None

    remainder = comment[len(class1_prefix):].strip()
    tokens = remainder.split()
    # Expected: tokens[0]="Frequency", tokens[1]="0.03659", tokens[2]="(African),"
    if len(tokens) < 3:
        return None, None

    frequency = tokens[1]
    population = tokens[2].strip('(),')
    return frequency, population


def create_brca_exchange_lookup(brca_exchange_output_tsv):
    """
    Read built_with_v4b.tsv and create a lookup dictionary.

    Returns a dictionary mapping variant keys (e.g., "BRCA2:c.-296C>T")
    to a dictionary containing V4_Popfreq, faf95_popmax_joint,
    faf95_popmax_population, faf95_popmax_allele_count, faf95_popmax_allele_number,
    enigma_class1_frequency, and enigma_class1_population values.
    """
    lookup = {}

    with open(brca_exchange_output_tsv, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        # Replace spaces with underscores in header


        for row in reader:
            gene_symbol = row["Gene_Symbol"]
            gnomADv4_id = row["Variant_id_GnomADv4"]
            pyhgvs_cdna = row["pyhgvs_cDNA"]
            v4_popfreq = row["Provisional_Evidence_Code_Gnomad_V4"]
            faf95_popmax_joint = row["faf95_popmax_joint_GnomADv4"]
            faf95_popmax_population = row["faf95_popmax_population_joint_GnomADv4"]
            comment = row["Comment_on_clinical_significance_ENIGMA"]

            # Get allele count and allele number values using the helper function
            faf95_popmax_allele_count, faf95_popmax_allele_number = get_allele_count_value(faf95_popmax_population, row)

            # Extract frequency and population from Class 1 comments
            enigma_class1_frequency, enigma_class1_population = parse_enigma_class1_comment(comment)

            # Skip rows with None or empty pyhgvs_cDNA
            if not pyhgvs_cdna or pyhgvs_cdna == "None" or pyhgvs_cdna == "-":
                continue

            # Create join key: replace part before first colon with gene_symbol
            if ':' in pyhgvs_cdna:
                # Split on first colon and reconstruct with gene_symbol
                parts = pyhgvs_cdna.split(':', 1)
                join_key = f"{gene_symbol}:{parts[1]}"
                lookup[join_key] = {
                    'gnomADv4_id': gnomADv4_id if gnomADv4_id else '',
                    'v4_popfreq': v4_popfreq if v4_popfreq else '',
                    'faf95_popmax_joint': faf95_popmax_joint if faf95_popmax_joint else '',
                    'faf95_popmax_population': faf95_popmax_population if faf95_popmax_population else '',
                    'faf95_popmax_allele_count': faf95_popmax_allele_count,
                    'faf95_popmax_allele_number': faf95_popmax_allele_number,
                    'enigma_class1_frequency': enigma_class1_frequency,
                    'enigma_class1_population': enigma_class1_population,
                }


    return lookup


def parse_arguments():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments with brca_exchange_output_tsv, enigma_tsv, and output_file.
    """
    parser = argparse.ArgumentParser(
        description='Merge V4_Popfreq data from built_with_v4b.tsv into ENIGMA_Johanna.tsv'
    )

    parser.add_argument(
        'brca_exchange_output_tsv',
        nargs='?',
        default='built_with_popfreq.tsv',
        help='Input file with V4_Popfreq data (default: built_with_popfreq.tsv)'
    )

    parser.add_argument(
        'enigma_tsv',
        nargs='?',
        default='ENIGMA_Johanna.xlsx',
        help='ENIGMA Excel spreadsheet to merge into (default: ENIGMA_Johanna.xlsx); reads from first sheet'
    )

    parser.add_argument(
        'output_file',
        nargs='?',
        default='ENIGMA_Johanna_merged.tsv',
        help='Output file path (default: ENIGMA_Johanna_merged.tsv)'
    )

    return parser.parse_args()


def merge_files(brca_exchange_output_tsv, enigma_tsv, output_file):
    """
    Merge V4_Popfreq, faf95_popmax_joint, faf95_popmax_population,
    faf95_popmax_allele_count, and faf95_popmax_allele_number data
    from brca_exchange_output_tsv into enigma_tsv.
    """
    csv.field_size_limit(sys.maxsize)
    # Create lookup dictionary from built_with_v4b.tsv
    print(f"Reading {brca_exchange_output_tsv}...")
    v4b_lookup = create_brca_exchange_lookup(brca_exchange_output_tsv)
    print(f"Created lookup with {len(v4b_lookup)} entries")

    # Read ENIGMA file and add columns
    print(f"Processing {enigma_tsv}...")
    matched = 0
    unmatched = 0

    if enigma_tsv.endswith('.xls'):
        xls_wb = xlrd.open_workbook(enigma_tsv)
        xls_ws = xls_wb.sheets()[0]
        rows_iter = (xls_ws.row_values(i) for i in range(xls_ws.nrows))
        wb = None
    else:
        wb = openpyxl.load_workbook(enigma_tsv, read_only=True, data_only=True)
        ws = wb.worksheets[0]
        rows_iter = ws.iter_rows(values_only=True)

    with open(output_file, 'w', encoding='utf-8', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')

        # Process header
        raw_header = next(rows_iter)
        header = [str(v) if v is not None else '' for v in raw_header[:13]]
        header = [field.rstrip().replace(' ', '_') for field in header]
        header.extend(['gnomADv4_id', 'faf95_popmax_joint_gnomADv4', 'faf95_popmax_population_joint_gnomADv4',
                      'faf95_popmax_allele_count_gnomADv4', 'faf95_popmax_allele_number_gnomADv4', 'V4_Popfreq',
                      'enigma_class1_frequency', 'enigma_class1_population'])
        writer.writerow(header)

        # Process data rows
        for raw_row in rows_iter:
            row = [str(v) if v is not None else '' for v in raw_row[:12]]
            if not any(row):
                continue

            row = [field.rstrip().replace(' ', '_') for field in row]

            variant = row[0]  # Column 1 (0-indexed as 0)

            # Look up values
            variant_data = v4b_lookup.get(variant, {})

            if variant_data:
                matched += 1
                gnomADv4_id = variant_data.get('gnomADv4_id', '')
                faf95_popmax_joint = variant_data.get('faf95_popmax_joint', '')
                faf95_popmax_population = variant_data.get('faf95_popmax_population', '')
                faf95_popmax_allele_count = variant_data.get('faf95_popmax_allele_count', '')
                faf95_popmax_allele_number = variant_data.get('faf95_popmax_allele_number', '')
                v4_popfreq = variant_data.get('v4_popfreq', '')
                enigma_class1_frequency = variant_data.get('enigma_class1_frequency') or ''
                enigma_class1_population = variant_data.get('enigma_class1_population') or ''
            else:
                unmatched += 1
                gnomADv4_id = ''
                faf95_popmax_joint = ''
                faf95_popmax_population = ''
                faf95_popmax_allele_count = ''
                faf95_popmax_allele_number = ''
                v4_popfreq = ''
                enigma_class1_frequency = ''
                enigma_class1_population = ''

            row.extend([gnomADv4_id, faf95_popmax_joint,
                        faf95_popmax_population, faf95_popmax_allele_count,
                        faf95_popmax_allele_number, v4_popfreq,
                        enigma_class1_frequency, enigma_class1_population])
            writer.writerow(row)

    if wb is not None:
        wb.close()
    print(f"Complete! Output written to {output_file}")
    print(f"Matched: {matched}, Unmatched: {unmatched}")

def main():
    args = parse_arguments()
    merge_files(args.brca_exchange_output_tsv, args.enigma_tsv, args.output_file)
    
if __name__ == '__main__':
    main()
