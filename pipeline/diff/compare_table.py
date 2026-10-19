#!/usr/bin/env python
"""
Compare a table across two PostgreSQL schemas.

Usage:
    python compare_table.py --table TABLE --old OLD_SCHEMA --new NEW_SCHEMA --output FILE
                            [--pk COL [COL ...]]

Stdout: counts of shared, deleted (in old not new), and added (in new not old) rows.

Output file: tab-separated rows for every column difference in shared rows:
    column_name  primary_key  old_value  new_value

Use --pk to override the primary key used for matching rows (useful when the
table's actual PK is a surrogate id that differs across schemas, e.g. report_gnomad).
Example: --pk VRS_Digest version data_type

Use --omit-columns to exclude columns from the diff (e.g. a surrogate id column
that legitimately differs across schemas and would otherwise show up as noise).
Example: --omit-columns id
"""

import argparse
import csv
import sys

import psycopg2
import psycopg2.extras


DB = dict(host='localhost', dbname='storage.pg', user='postgres', password='postgres')


def _resolve_base_table(cur, schema, table):
    """If table is a view, return the underlying base table name; else return table."""
    cur.execute("""
        SELECT view_definition
        FROM information_schema.views
        WHERE table_schema = %s AND table_name = %s
    """, [schema, table])
    row = cur.fetchone()
    if row is None:
        return table
    # Extract base table name from "SELECT ... FROM schema.tablename"
    import re
    m = re.search(r'FROM\s+"?(\w+)"?\."?(\w+)"?', row['view_definition'], re.IGNORECASE)
    return m.group(2) if m else table


def get_pk_columns(cur, schema, table):
    base = _resolve_base_table(cur, schema, table)
    cur.execute("""
        SELECT kcu.column_name
        FROM information_schema.table_constraints tc
        JOIN information_schema.key_column_usage kcu
            ON tc.constraint_name = kcu.constraint_name
           AND tc.table_schema    = kcu.table_schema
        WHERE tc.constraint_type = 'PRIMARY KEY'
          AND tc.table_schema    = %s
          AND tc.table_name      = %s
        ORDER BY kcu.ordinal_position
    """, [schema, base])
    cols = [row['column_name'] for row in cur.fetchall()]
    if not cols:
        sys.exit(f"No primary key found for {schema}.{table}")
    return cols


def get_columns(cur, schema, table):
    cur.execute("""
        SELECT column_name
        FROM information_schema.columns
        WHERE table_schema = %s AND table_name = %s
        ORDER BY ordinal_position
    """, [schema, table])
    return [row['column_name'] for row in cur.fetchall()]


def pk_expr(pk_cols, alias):
    """Return SQL fragment: alias.col1, alias.col2, ..."""
    return ', '.join(f'{alias}."{c}"' for c in pk_cols)


def join_cond(pk_cols, alias_old, alias_new):
    return ' AND '.join(
        f'{alias_old}."{c}" = {alias_new}."{c}"' for c in pk_cols
    )


def pk_value(row, pk_cols):
    vals = [str(row[c]) for c in pk_cols]
    return '|'.join(vals) if len(vals) > 1 else vals[0]


def fetch_rows(cur, schema, table):
    cur.execute(f'SELECT * FROM "{schema}"."{table}"')
    return {pk_value(r, cur.description): r for r in cur.fetchall()}  # placeholder


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--table',  required=True, help='Table name')
    parser.add_argument('--old',    required=True, help='Old schema name')
    parser.add_argument('--new',    required=True, help='New schema name')
    parser.add_argument('--output', required=True, help='Detailed diff output file')
    parser.add_argument('--pk', nargs='+', metavar='COL',
                        help='Override primary key columns used for row matching')
    parser.add_argument('--omit-columns', nargs='+', metavar='COL', default=[],
                        help='Columns to exclude from the diff')
    args = parser.parse_args()

    conn = psycopg2.connect(**DB)
    conn.autocommit = True

    with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
        pk_cols = args.pk if args.pk else get_pk_columns(cur, args.old, args.table)
        all_cols = get_columns(cur, args.old, args.table)
        data_cols = [c for c in all_cols if c not in pk_cols and c not in args.omit_columns]

        pk_join = join_cond(pk_cols, 'o', 'n')
        old_pk_null = ' AND '.join(f'o."{c}" IS NULL' for c in pk_cols)
        new_pk_null = ' AND '.join(f'n."{c}" IS NULL' for c in pk_cols)

        # --- counts ---
        cur.execute(f"""
            SELECT COUNT(*) FROM "{args.old}"."{args.table}" o
            JOIN "{args.new}"."{args.table}" n ON {pk_join}
        """)
        n_shared = cur.fetchone()['count']

        cur.execute(f"""
            SELECT COUNT(*) FROM "{args.old}"."{args.table}" o
            LEFT JOIN "{args.new}"."{args.table}" n ON {pk_join}
            WHERE {new_pk_null}
        """)
        n_deleted = cur.fetchone()['count']

        cur.execute(f"""
            SELECT COUNT(*) FROM "{args.new}"."{args.table}" n
            LEFT JOIN "{args.old}"."{args.table}" o ON {pk_join}
            WHERE {old_pk_null}
        """)
        n_added = cur.fetchone()['count']

        print(f"shared:  {n_shared}")
        print(f"deleted: {n_deleted}  (in {args.old}, not in {args.new})")
        print(f"added:   {n_added}  (in {args.new}, not in {args.old})")

        # --- detailed diff of shared rows ---
        # Fetch both sides of shared rows in one query
        old_sel = ', '.join(f'o."{c}" AS "old_{c}"' for c in data_cols)
        new_sel = ', '.join(f'n."{c}" AS "new_{c}"' for c in data_cols)
        pk_sel  = ', '.join(f'o."{c}" AS "{c}"' for c in pk_cols)

        cur.execute(f"""
            SELECT {pk_sel}, {old_sel}, {new_sel}
            FROM "{args.old}"."{args.table}" o
            JOIN "{args.new}"."{args.table}" n ON {pk_join}
        """)

        with open(args.output, 'w', newline='') as fh:
            writer = csv.writer(fh, delimiter='\t')
            writer.writerow(['column', 'primary_key', 'old_value', 'new_value'])
            n_diff_rows = 0
            for row in cur:
                pk = pk_value(dict(row), pk_cols)
                for col in data_cols:
                    old_val = row[f'old_{col}']
                    new_val = row[f'new_{col}']
                    if old_val != new_val:
                        writer.writerow([col, pk, old_val, new_val])
                        n_diff_rows += 1

        print(f"differences in shared rows written to {args.output} ({n_diff_rows} differing values)")

    conn.close()


if __name__ == '__main__':
    main()
