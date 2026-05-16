"""
Populate enigma_domain with the ENIGMA Consortium functional domains of
potential clinical importance for BRCA1 and BRCA2.

The table is truncated and reloaded on each run so the data stays in sync
with the authoritative source (website/js/SplicingData.js).
"""

import click
import psycopg2

DOMAINS = [
    # (gene, name, chrom, assembly, start, end)
    ('BRCA1', 'ring',          '17', 'GRCh38', 43104260, 43124093),
    ('BRCA1', 'coiled-coil',   '17', 'GRCh38', 43082489, 43090958),
    ('BRCA1', 'brct',          '17', 'GRCh38', 43045699, 43070966),
    ('BRCA2', 'palb2 binding', '13', 'GRCh38', 32316488, 32319129),
    ('BRCA2', 'dna binding',   '13', 'GRCh38', 32356433, 32396954),
]


@click.command()
@click.option('--db-url', default='postgresql://postgres:postgres@localhost/storage.pg',
              envvar='PIPELINE_DB_URL', show_default=True,
              help='PostgreSQL connection URL for the pipeline DB')
@click.option('--schema', default='pipeline', show_default=True,
              help='Schema containing enigma_domain')
def main(db_url, schema):
    conn = psycopg2.connect(db_url, options=f'-c search_path={schema}')
    try:
        with conn.cursor() as cur:
            cur.execute('TRUNCATE enigma_domain RESTART IDENTITY')
            cur.executemany(
                """INSERT INTO enigma_domain (gene, name, chrom, assembly, start, "end")
                   VALUES (%s, %s, %s, %s, %s, %s)""",
                DOMAINS,
            )
        conn.commit()
        print(f'Loaded {len(DOMAINS)} ENIGMA domain rows.')
    finally:
        conn.close()


if __name__ == '__main__':
    main()
