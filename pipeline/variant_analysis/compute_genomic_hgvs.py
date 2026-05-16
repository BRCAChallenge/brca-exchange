"""
Populate genomic_coordinates.hgvs with normalized genomic HGVS strings.

For every row the genomic HGVS string is derived from the chr/pos/ref/alt
fields using the biocommons hgvs library and the GRCh38 contig accession map.
Normalization is attempted; if it fails the unnormalized form is stored.

Requires:
  - UTA_DB_URL env var pointing to a local UTA PostgreSQL instance
  - HGVS_SEQREPO_DIR env var pointing to a local seqrepo directory
"""

import logging
import os

import click
import psycopg2

from common.hgvs_utils import HgvsWrapper
from common.variant_utils import VCFVariant

BATCH_SIZE = 500


def _genomic_hgvs(wrapper, contig_map, chrom, pos, ref, alt):
    chrom = str(chrom).removeprefix('chr')
    vcf_var = VCFVariant(chrom, int(pos), ref, alt)
    hgvs_obj = vcf_var.to_hgvs_obj(contig_map)
    try:
        return str(wrapper.hgvs_norm.normalize(hgvs_obj))
    except Exception:
        return str(hgvs_obj)


@click.command()
@click.option('--db-url', default='postgresql://postgres:postgres@localhost/storage.pg',
              envvar='PIPELINE_DB_URL', show_default=True,
              help='PostgreSQL connection URL for the pipeline DB')
@click.option('--schema', default='pipeline', show_default=True,
              help='Schema containing genomic_coordinates')
def main(db_url, schema):
    logging.basicConfig(level=logging.WARNING,
                        format='%(levelname)s %(message)s')

    wrapper = HgvsWrapper()
    contig_map = wrapper.contig_maps[HgvsWrapper.GRCh38_Assem]

    conn = psycopg2.connect(db_url, options=f'-c search_path={schema}')
    try:
        with conn.cursor() as cur:
            cur.execute('SELECT id, chr, pos, ref, alt FROM genomic_coordinates')
            rows = cur.fetchall()

        print(f'Computing genomic HGVS for {len(rows)} rows ...')

        updates = []
        failed = 0
        for row_id, chrom, pos, ref, alt in rows:
            try:
                hgvs_str = _genomic_hgvs(wrapper, contig_map, chrom, pos, ref, alt)
                updates.append((hgvs_str, row_id))
            except Exception as e:
                logging.warning('HGVS failed for id=%s %s:%s:%s>%s: %s',
                                row_id, chrom, pos, ref, alt, e)
                failed += 1

        with conn.cursor() as cur:
            for i in range(0, len(updates), BATCH_SIZE):
                batch = updates[i : i + BATCH_SIZE]
                cur.executemany(
                    'UPDATE genomic_coordinates SET hgvs = %s WHERE id = %s',
                    batch,
                )
        conn.commit()
        print(f'Updated {len(updates)} rows ({failed} failures).')
    finally:
        conn.close()


if __name__ == '__main__':
    main()
