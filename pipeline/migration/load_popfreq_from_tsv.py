"""
Load Provisional_Evidence_Code_Popfreq / _Description_Popfreq from a release TSV
into pipeline.analysis_provisional_evidence_codes.

Usage:
    python load_popfreq_from_tsv.py --tsv built_with_change_types.tsv

Columns used (1-based, matching the built_with_change_types.tsv layout):
    455  VR_ID                               → VRS_Digest
    466  Provisional_Evidence_Code_Popfreq   → popfreq_code
    467  Provisional_Evidence_Description_Popfreq → popfreq_description

method_name is fixed to 'popfreq_1.1'.
Rows whose VR_ID is empty are skipped.
Existing rows (same VRS_Digest + method_name) are updated in place.
"""

import argparse
import csv
import sys

import psycopg2
import psycopg2.extras

DB_URL = 'postgresql://postgres:postgres@localhost/storage.pg'
SCHEMA = 'pipeline'
METHOD_NAME = 'popfreq_1.1'
BATCH = 2000

COL_VRS = 'VR_ID'
COL_CODE = 'Provisional_Evidence_Code_Popfreq'
COL_DESC = 'Provisional_Evidence_Description_Popfreq'

csv.field_size_limit(10 ** 7)


def load(tsv_path, db_url, dry_run):
    rows = []
    skipped = 0
    with open(tsv_path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for rec in reader:
            vrs = rec.get(COL_VRS, '').strip()
            if not vrs:
                skipped += 1
                continue
            rows.append((
                vrs,
                rec.get(COL_CODE, '').strip() or None,
                rec.get(COL_DESC, '').strip() or None,
                METHOD_NAME,
            ))

    print(f'Read {len(rows)} rows ({skipped} skipped — empty VR_ID)')

    if dry_run:
        print('Dry run — first 3 rows:')
        for r in rows[:3]:
            print(' ', r)
        return

    conn = psycopg2.connect(db_url, options=f'-c search_path={SCHEMA}')
    try:
        with conn.cursor() as cur:
            cur.execute('SELECT "VRS_Digest" FROM variant')
            valid_digests = {r[0] for r in cur.fetchall()}
        before = len(rows)
        rows = [r for r in rows if r[0] in valid_digests]
        dropped = before - len(rows)
        if dropped:
            print(f'Dropped {dropped} rows whose VRS_Digest is not in the variant table')

        with conn.cursor() as cur:
            for i in range(0, len(rows), BATCH):
                psycopg2.extras.execute_values(
                    cur,
                    """INSERT INTO analysis_provisional_evidence_codes
                           ("VRS_Digest", popfreq_code, popfreq_description, method_name)
                       VALUES %s
                       ON CONFLICT ("VRS_Digest", method_name) DO UPDATE
                           SET popfreq_code        = EXCLUDED.popfreq_code,
                               popfreq_description = EXCLUDED.popfreq_description""",
                    rows[i:i + BATCH],
                )
                if (i // BATCH) % 10 == 0:
                    print(f'  {min(i + BATCH, len(rows))}/{len(rows)} rows written…')
        conn.commit()
        print(f'Done. Inserted/updated {len(rows)} rows with method_name={METHOD_NAME!r}.')
    finally:
        conn.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--tsv', required=True, help='Path to built_with_change_types.tsv')
    parser.add_argument('--db-url', default=DB_URL, help='PostgreSQL connection URL')
    parser.add_argument('--dry-run', action='store_true',
                        help='Parse TSV and print sample output without writing to DB')
    args = parser.parse_args()
    load(args.tsv, args.db_url, args.dry_run)


if __name__ == '__main__':
    main()
