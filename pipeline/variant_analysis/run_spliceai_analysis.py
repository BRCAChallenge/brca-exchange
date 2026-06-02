"""
Populate analysis_spliceai from a SpliceAI-annotated VCF.

Parses DS_AG, DS_AL, DS_DG, DS_DL, DP_AG, DP_AL, DP_DG, DP_DL from the
SpliceAI INFO field and joins to the database by GRCh38 chr/pos/ref/alt.
result = max(DS_AG, DS_AL, DS_DG, DS_DL).
Variants already in analysis_spliceai are skipped unless --overwrite is set.
"""

import logging

import click
import psycopg2
import psycopg2.extras

DB_BATCH = 500

log = logging.getLogger(__name__)

# SpliceAI INFO format: SpliceAI=ALLELE|GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
_SPLICEAI_FIELDS = ['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL',
                    'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL']


def _parse_spliceai_vcf(path):
    """Return dict mapping (chr, pos, ref, alt) -> score dict."""
    scores = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 8:
                continue
            chrom, pos, ref, alt, info = cols[0], cols[1], cols[3], cols[4], cols[7]
            for field in info.split(';'):
                if not field.startswith('SpliceAI='):
                    continue
                parts = field[len('SpliceAI='):].split('|')
                if len(parts) < 10:
                    continue
                # parts: ALLELE GENE DS_AG DS_AL DS_DG DS_DL DP_AG DP_AL DP_DG DP_DL
                ds_values = []
                row = {}
                for i, name in enumerate(_SPLICEAI_FIELDS):
                    val = parts[i + 2]  # skip ALLELE and GENE
                    row[name] = val
                    if name.startswith('DS_'):
                        try:
                            ds_values.append(float(val))
                        except ValueError:
                            pass
                row['result'] = str(max(ds_values)) if ds_values else None
                scores[(chrom, pos, ref, alt)] = row
                break
    return scores


@click.command()
@click.option('--db-url', default='postgresql://postgres:postgres@localhost/storage.pg',
              envvar='PIPELINE_DB_URL', show_default=True)
@click.option('--schema', default='pipeline', show_default=True)
@click.option('--spliceai-vcf', required=True, type=click.Path(exists=True),
              help='Path to the SpliceAI-annotated VCF')
@click.option('--overwrite', is_flag=True, default=False,
              help='Re-score variants already in analysis_spliceai')
def main(db_url, schema, spliceai_vcf, overwrite):
    logging.basicConfig(level=logging.WARNING, format='%(levelname)s %(message)s')

    print(f'Parsing SpliceAI VCF: {spliceai_vcf} ...')
    scores = _parse_spliceai_vcf(spliceai_vcf)
    print(f'  {len(scores)} SpliceAI records loaded from VCF.')

    conn = psycopg2.connect(db_url, options=f'-c search_path={schema}')
    try:
        with conn.cursor() as cur:
            if overwrite:
                cur.execute("""
                    SELECT gc."VRS_Digest_id", gc.chr, gc.pos, gc.ref, gc.alt
                    FROM genomic_coordinates gc
                    WHERE gc.assembly = 'GRCh38'
                """)
            else:
                cur.execute("""
                    SELECT gc."VRS_Digest_id", gc.chr, gc.pos, gc.ref, gc.alt
                    FROM genomic_coordinates gc
                    LEFT JOIN analysis_spliceai sa ON sa."VRS_Digest" = gc."VRS_Digest_id"
                    WHERE gc.assembly = 'GRCh38'
                      AND sa."VRS_Digest" IS NULL
                """)
            coord_rows = cur.fetchall()

        print(f'Matching {len(coord_rows)} variants against VCF scores ...')
        inserts = []
        for vrs_digest, chr_, pos, ref, alt in coord_rows:
            chr_bare = chr_[3:] if chr_.startswith('chr') else chr_
            row = scores.get((chr_bare, pos, ref, alt))
            if row is not None:
                inserts.append((
                    vrs_digest,
                    row['DS_AG'], row['DS_AL'], row['DS_DG'], row['DS_DL'],
                    row['DP_AG'], row['DP_AL'], row['DP_DG'], row['DP_DL'],
                    row['result'],
                ))

        print(f'Writing {len(inserts)} rows to analysis_spliceai ...')
        with conn.cursor() as cur:
            for i in range(0, len(inserts), DB_BATCH):
                psycopg2.extras.execute_values(
                    cur,
                    """INSERT INTO analysis_spliceai
                           ("VRS_Digest", "DS_AG", "DS_AL", "DS_DG", "DS_DL",
                            "DP_AG", "DP_AL", "DP_DG", "DP_DL", result)
                       VALUES %s
                       ON CONFLICT ("VRS_Digest") DO UPDATE SET
                           "DS_AG" = EXCLUDED."DS_AG",
                           "DS_AL" = EXCLUDED."DS_AL",
                           "DS_DG" = EXCLUDED."DS_DG",
                           "DS_DL" = EXCLUDED."DS_DL",
                           "DP_AG" = EXCLUDED."DP_AG",
                           "DP_AL" = EXCLUDED."DP_AL",
                           "DP_DG" = EXCLUDED."DP_DG",
                           "DP_DL" = EXCLUDED."DP_DL",
                           result  = EXCLUDED.result""",
                    inserts[i : i + DB_BATCH],
                )
        conn.commit()
        print(f'Done. Inserted/updated {len(inserts)} rows.')
    finally:
        conn.close()


if __name__ == '__main__':
    main()
