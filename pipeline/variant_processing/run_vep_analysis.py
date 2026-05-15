"""
Run VEP on all variants and populate analysis_vep.

For each variant with a GRCh38 genomic HGVS in genomic_coordinates,
queries the local VEP server in batches and stores variant_class in
the analysis_vep table.

Variants that already have a row in analysis_vep are skipped unless
--overwrite is set.
"""

import logging

import click
import psycopg2
import psycopg2.extras
import requests

BATCH_SIZE   = 200   # variants per VEP call
DB_BATCH     = 500   # rows per INSERT batch

log = logging.getLogger(__name__)


def _query_vep_batch(session, vep_url, hgvs_list):
    """POST a batch of HGVS strings to the VEP server.
    Returns dict {hgvs: [records]} or {hgvs: {'error': ...}}.
    """
    resp = session.post(f'{vep_url}/vep/hgvs/batch',
                        json=hgvs_list, timeout=600)
    resp.raise_for_status()
    return resp.json()


@click.command()
@click.option('--db-url', default='postgresql://postgres:postgres@localhost/storage.pg',
              envvar='PIPELINE_DB_URL', show_default=True)
@click.option('--schema', default='pipeline', show_default=True)
@click.option('--vep-url', default='http://localhost:8888', show_default=True,
              envvar='VEP_SERVER_URL', help='Base URL of the local VEP REST server')
@click.option('--overwrite', is_flag=True, default=False,
              help='Re-annotate variants that already have a row in analysis_vep')
def main(db_url, schema, vep_url, overwrite):
    logging.basicConfig(level=logging.WARNING, format='%(levelname)s %(message)s')

    conn = psycopg2.connect(db_url, options=f'-c search_path={schema}')
    try:
        with conn.cursor() as cur:
            if overwrite:
                cur.execute("""
                    SELECT gc."VRS_Digest_id", gc.hgvs
                    FROM genomic_coordinates gc
                    WHERE gc.assembly = 'GRCh38'
                      AND gc.hgvs IS NOT NULL AND gc.hgvs <> ''
                """)
            else:
                cur.execute("""
                    SELECT gc."VRS_Digest_id", gc.hgvs
                    FROM genomic_coordinates gc
                    LEFT JOIN analysis_vep av ON av."VRS_Digest" = gc."VRS_Digest_id"
                    WHERE gc.assembly = 'GRCh38'
                      AND gc.hgvs IS NOT NULL AND gc.hgvs <> ''
                      AND av."VRS_Digest" IS NULL
                """)
            rows = cur.fetchall()

        total = len(rows)
        print(f'Running VEP on {total} variants in batches of {BATCH_SIZE} ...')

        session = requests.Session()
        inserts = []
        errors  = 0

        for batch_start in range(0, total, BATCH_SIZE):
            batch = rows[batch_start : batch_start + BATCH_SIZE]
            digests  = [r[0] for r in batch]
            hgvs_map = {r[1]: r[0] for r in batch}   # hgvs → VRS_Digest

            if (batch_start // BATCH_SIZE) % 10 == 0 or batch_start + BATCH_SIZE >= total:
                done = min(batch_start + BATCH_SIZE, total)
                print(f'  {done}/{total}  ({len(inserts)} queued, {errors} errors)')

            try:
                results = _query_vep_batch(session, vep_url, list(hgvs_map.keys()))
            except Exception as e:
                log.warning('Batch starting at %d failed: %s', batch_start, e)
                errors += len(batch)
                continue

            for hgvs, records in results.items():
                vrs = hgvs_map.get(hgvs)
                if not vrs:
                    continue
                if isinstance(records, dict) and 'error' in records:
                    log.warning('VEP error for %s: %s', hgvs, records['error'])
                    errors += 1
                    continue
                variant_class = records[0].get('variant_class') if records else None
                inserts.append((vrs, variant_class))

        print(f'Writing {len(inserts)} rows to analysis_vep ...')
        with conn.cursor() as cur:
            for i in range(0, len(inserts), DB_BATCH):
                psycopg2.extras.execute_values(
                    cur,
                    """INSERT INTO analysis_vep ("VRS_Digest", variant_class)
                       VALUES %s
                       ON CONFLICT ("VRS_Digest") DO UPDATE
                         SET variant_class = EXCLUDED.variant_class""",
                    inserts[i : i + DB_BATCH],
                )
        conn.commit()
        print(f'Done. Inserted/updated {len(inserts)}, errors {errors}.')
    finally:
        conn.close()


if __name__ == '__main__':
    main()
