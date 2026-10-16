#!/usr/bin/env python
"""
Compare tables across two schemas.

Usage:
    python compare_all_tables.py --old OLD_SCHEMA --new NEW_SCHEMA --output-dir DIR

Runs compare_table.py for each table and writes output files to DIR:
    variant.detail.txt
    variant_genomic_coordinates.detail.txt
    variant_gnomad.detail.txt
    report_gnomad.detail.txt
    variant_lovd.detail.txt
    report_lovd.detail.txt
    variant_clinvar.detail.txt
    report_clinvar.detail.txt
    variant_enigma.detail.txt
    variant_exlovd.detail.txt
    variant_other.detail.txt
    analysis_bayesdel.detail.txt
    analysis_priors.detail.txt
    analysis_provisional_evidence_codes.detail.txt
    summary.txt

variant uses VRS_Digest as its natural primary key.

variant_genomic_coordinates uses (VRS_Digest, assembly) as the natural key for
matching rows, since its surrogate id column differs across schemas. That same
id column is also omitted from the diff itself, since it always differs.

report_gnomad uses (VRS_Digest, version, data_type) as the natural key for
matching rows, since its surrogate id column differs across schemas. That same
id column is also omitted from the diff itself, since it always differs.

report_lovd uses (VRS_Digest, Submission_ID) as the natural key, for the
same reason (surrogate id column differs across schemas); that id column is
likewise omitted from the diff.

report_clinvar uses (VRS_Digest, SCV) as the natural key, for the same
reason; SCV (the ClinVar submission accession) is unique on its own, but
VRS_Digest is included for consistency with the other report_* tables.
That surrogate id column is likewise omitted from the diff.

variant_enigma (despite the name, a report_* table with a surrogate id: one
ENIGMA row per variant in practice, via a plain ForeignKey rather than
OneToOne) uses VRS_Digest alone as the natural key, since it's unique on
its own; the id column is likewise omitted from the diff.

summary.txt collects, for each table, the stdout (shared/deleted/added counts
and the differing-values line) of the compare_table.py run for that table.
"""

import argparse
import os
import subprocess
import sys

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_COMPARE = os.path.join(_SCRIPT_DIR, 'compare_table.py')

TABLES = [
    dict(table='variant',              pk=None, omit=['Synonyms', 'VRS', 'ensembl_cdna', 'ensembl_protein', 'title']),
    dict(table='variant_genomic_coordinates',  pk=['VRS_Digest', 'assembly'], omit=['id']),
    dict(table='variant_gnomad', pk=None, omit=None),
    dict(table='report_gnomad',  pk=['VRS_Digest', 'version', 'data_type'], omit=['id']),
    dict(table='variant_lovd',   pk=None, omit=None),
    dict(table='report_lovd',    pk=['VRS_Digest', 'Submission_ID'], omit=['id']),
    dict(table='variant_clinvar', pk=None, omit=None),
    dict(table='report_clinvar',  pk=['VRS_Digest', 'SCV'], omit=['id']),
    dict(table='variant_enigma',  pk=['VRS_Digest'], omit=['id']),
    dict(table='variant_exlovd',  pk=None, omit=None),
    dict(table='variant_other',        pk=None, omit=None),
    dict(table='analysis_bayesdel',   pk=None, omit=None),
    dict(table='analysis_priors',     pk=None, omit=None),
    dict(table='analysis_provisional_evidence_codes', pk=['VRS_Digest', 'method_name'], omit=['id']),
]


def run(old, new, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    python = sys.executable
    summary_path = os.path.join(output_dir, 'summary.txt')
    with open(summary_path, 'w') as summary:
        for spec in TABLES:
            table = spec['table']
            output = os.path.join(output_dir, f'{table}.detail.txt')
            cmd = [python, _COMPARE,
                   '--table', table,
                   '--old', old,
                   '--new', new,
                   '--output', output]
            if spec['pk']:
                cmd += ['--pk'] + spec['pk']
            if spec['omit']:
                cmd += ['--omit-columns'] + spec['omit']
            print(f'\n=== {table} ===')
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(result.stdout, end='')
            summary.write(f'=== {table} ===\n{result.stdout}\n')
    print(f'\nSummary written to {summary_path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--old',        required=True, help='Old schema name')
    parser.add_argument('--new',        required=True, help='New schema name')
    parser.add_argument('--output-dir', required=True, help='Directory for output files')
    args = parser.parse_args()
    run(args.old, args.new, args.output_dir)


if __name__ == '__main__':
    main()
