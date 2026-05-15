"""
Populate analysis_bayesdel from a Victor-annotated VCF.

Reads BayesDel_nsfp33a_noAF from the VCF INFO field and joins to the
database by GRCh38 chr/pos/ref/alt to obtain VRS_Digest.
Variants already in analysis_bayesdel are skipped unless --overwrite is set.
"""

import logging

import click
import psycopg2
import psycopg2.extras

DB_BATCH = 500

log = logging.getLogger(__name__)


def _parse_victor_vcf(path):
    """Return dict mapping (chr, pos, ref, alt) -> BayesDel_nsfp33a_noAF string."""
    scores = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 8:
                continue
            chrom, pos, ref, alt, info = cols[0], cols[1], cols[3], cols[4], cols[7]
            score = None
            for field in info.split(';'):
                if field.startswith('BayesDel_nsfp33a_noAF='):
                    score = field.split('=', 1)[1]
                    break
            if score is not None and score not in ('.', ''):
                scores[(chrom, pos, ref, alt)] = score
    return scores


@click.command()
@click.option('--db-url', default='postgresql://postgres:postgres@localhost/storage.pg',
              envvar='PIPELINE_DB_URL', show_default=True)
@click.option('--schema', default='pipeline', show_default=True)
@click.option('--victor-vcf', required=True, type=click.Path(exists=True),
              help='Path to the Victor-annotated VCF (BRCA.ann.all.vcf)')
@click.option('--overwrite', is_flag=True, default=False,
              help='Re-score variants already in analysis_bayesdel')
def main(db_url, schema, victor_vcf, overwrite):
    logging.basicConfig(level=logging.WARNING, format='%(levelname)s %(message)s')

    print(f'Parsing Victor VCF: {victor_vcf} ...')
    scores = _parse_victor_vcf(victor_vcf)
    print(f'  {len(scores)} BayesDel scores loaded from VCF.')

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
                    LEFT JOIN analysis_bayesdel ab ON ab."VRS_Digest" = gc."VRS_Digest_id"
                    WHERE gc.assembly = 'GRCh38'
                      AND ab."VRS_Digest" IS NULL
                """)
            coord_rows = cur.fetchall()

        print(f'Matching {len(coord_rows)} variants against VCF scores ...')
        inserts = []
        for vrs_digest, chr_, pos, ref, alt in coord_rows:
            chr_bare = chr_[3:] if chr_.startswith('chr') else chr_
            score = scores.get((chr_bare, pos, ref, alt))
            if score is not None:
                inserts.append((vrs_digest, score))

        print(f'Writing {len(inserts)} rows to analysis_bayesdel ...')
        with conn.cursor() as cur:
            for i in range(0, len(inserts), DB_BATCH):
                psycopg2.extras.execute_values(
                    cur,
                    """INSERT INTO analysis_bayesdel ("VRS_Digest", "BayesDel_nsfp33a_noAF")
                       VALUES %s
                       ON CONFLICT ("VRS_Digest") DO UPDATE
                         SET "BayesDel_nsfp33a_noAF" = EXCLUDED."BayesDel_nsfp33a_noAF" """,
                    inserts[i : i + DB_BATCH],
                )
        conn.commit()
        print(f'Done. Inserted/updated {len(inserts)} rows.')
    finally:
        conn.close()


if __name__ == '__main__':
    main()
