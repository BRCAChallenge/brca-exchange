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


def create_brca_exchange_lookup(brca_exchange_output_tsv):
    """
    Read built_with_v4b.tsv and create a lookup dictionary.

    Returns a dictionary mapping variant keys (e.g., "BRCA2:c.-296C>T")
    to a dictionary containing V4_Popfreq, faf95_popmax_joint,
    faf95_popmax_population, faf95_popmax_allele_count, and faf95_popmax_allele_number values.
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

            # Get allele count and allele number values using the helper function
            faf95_popmax_allele_count, faf95_popmax_allele_number = get_allele_count_value(faf95_popmax_population, row)

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
                    'faf95_popmax_allele_number': faf95_popmax_allele_number
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
        default='built_with_v4_verify.tsv',
        help='Input file with V4_Popfreq data (default: built_with_v4b.tsv)'
    )

    parser.add_argument(
        'enigma_tsv',
        nargs='?',
        default='ENIGMA_Johanna.tsv',
        help='ENIGMA file to merge into (default: ENIGMA_Johanna.tsv)'
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

    with open(enigma_tsv, 'r', encoding='utf-8') as infile, \
         open(output_file, 'w', encoding='utf-8', newline='') as outfile:

        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        # Process header
        header = next(reader)
        # Keep only the first 12 columns
        header = header[:12]
        # Trim trailing whitespace from column 12 (index 11) and replace spaces with underscores
        header = [field.rstrip().replace(' ', '_') if i >= 11 else field.replace(' ', '_')
                  for i, field in enumerate(header)]
        header.extend(['gnomADv4_id', 'faf95_popmax_joint_gnomADv4', 'faf95_popmax_population_joint_gnomADv4',
                      'faf95_popmax_allele_count_gnomADv4', 'faf95_popmax_allele_number_gnomADv4', 'V4_Popfreq'])
        writer.writerow(header)

        # Process data rows
        for row in reader:
            if len(row) == 0:
                continue

            # Keep only the first 11 columns, or pad to 11 columns if shorter
            row = row[:11]
            # Pad with empty strings if row has fewer than 11 values
            while len(row) < 11:
                row.append('')

            # Trim trailing whitespace from column 12 (index 11) and replace spaces with underscores
            row = [field.rstrip().replace(' ', '_') if i >= 10 else field.replace(' ', '_')
                   for i, field in enumerate(row)]

            variant = row[0]  # Column 1 (0-indexed as 0)

            # Look up values
            variant_data = v4b_lookup.get(variant, {})

            if variant_data:
                matched += 1
                gnomADv4_id = variant_data.get('gnomADv4_id', '')
                print("Setting output variable to ", gnomADv4_id)
                faf95_popmax_joint = variant_data.get('faf95_popmax_joint', '')
                faf95_popmax_population = variant_data.get('faf95_popmax_population', '')
                faf95_popmax_allele_count = variant_data.get('faf95_popmax_allele_count', '')
                faf95_popmax_allele_number = variant_data.get('faf95_popmax_allele_number', '')
                v4_popfreq = variant_data.get('v4_popfreq', '')
            else:
                unmatched += 1
                faf95_popmax_joint = ''
                faf95_popmax_population = ''
                faf95_popmax_allele_count = ''
                faf95_popmax_allele_number = ''
                v4_popfreq = ''

            row.extend([gnomADv4_id, faf95_popmax_joint,
                        faf95_popmax_population, faf95_popmax_allele_count,
                        faf95_popmax_allele_number, v4_popfreq])
            writer.writerow(row)

    print(f"Complete! Output written to {output_file}")
    print(f"Matched: {matched}, Unmatched: {unmatched}")

def main():
    args = parse_arguments()
    merge_files(args.brca_exchange_output_tsv, args.enigma_tsv, args.output_file)
    
if __name__ == '__main__':
    main()
