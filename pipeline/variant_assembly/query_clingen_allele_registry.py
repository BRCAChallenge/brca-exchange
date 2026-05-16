"""
Query the ClinGen Allele Registry for each variant using the GRCh38 genomic
HGVS string, then populate variant fields:

  - CA_ID           ← suffix of resp['@id']             e.g. "CA276568"
  - title           ← resp['communityStandardTitle'][0]
  - HGVS_cDNA       ← MANE Select nucleotide RefSeq hgvs
  - ensembl_cdna    ← MANE Select nucleotide Ensembl hgvs
  - HGVS_Protein    ← MANE Select protein RefSeq hgvs
  - ensembl_protein ← MANE Select protein Ensembl hgvs
  - Synonyms        ← '|'-joined hgvs + protein hgvs from all non-MANE transcript alleles

Also inserts a GRCh37 row into genomic_coordinates for each genomicAllele
whose referenceGenome is 'GRCh37' (skipped if row already exists).

Variants that already have a CA_ID are skipped unless --overwrite is set.
"""

import logging
import time

import click
import psycopg2
import psycopg2.extras
import requests

CLINGEN_URL = 'https://reg.clinicalgenome.org/allele'
BATCH_SIZE = 500
REQUEST_DELAY = 0.1   # seconds between requests
MAX_RETRIES = 3
RETRY_BACKOFF = 2.0   # seconds; doubles on each retry


def _query(session, hgvs):
    """Return parsed JSON from ClinGen, or None on 404, or raise."""
    delay = RETRY_BACKOFF
    for attempt in range(MAX_RETRIES):
        try:
            resp = session.get(CLINGEN_URL, params={'hgvs': hgvs}, timeout=15)
            if resp.status_code == 200:
                return resp.json()
            if resp.status_code == 429:
                time.sleep(delay)
                delay *= 2
                continue
            if resp.status_code == 404:
                return None
            resp.raise_for_status()
        except requests.RequestException as e:
            if attempt == MAX_RETRIES - 1:
                raise
            logging.warning('Retry %d for %s: %s', attempt + 1, hgvs, e)
            time.sleep(delay)
            delay *= 2
    return None


def _mane_select(data):
    """Return (nuc_refseq, nuc_ensembl, prot_refseq, prot_ensembl) from the
    MANE Select transcript allele, or (None, None, None, None)."""
    for ta in data.get('transcriptAlleles', []):
        mane = ta.get('MANE')
        if not mane or mane.get('maneStatus') != 'MANE Select':
            continue
        nuc  = mane.get('nucleotide', {})
        prot = mane.get('protein',    {})
        return (
            nuc.get('RefSeq',  {}).get('hgvs'),
            nuc.get('Ensembl', {}).get('hgvs'),
            prot.get('RefSeq',  {}).get('hgvs'),
            prot.get('Ensembl', {}).get('hgvs'),
        )
    return None, None, None, None


def _synonyms(data):
    """Return a '|'-joined string of all hgvs and protein hgvs from transcript
    alleles that do not have a MANE object."""
    terms = []
    for ta in data.get('transcriptAlleles', []):
        if 'MANE' in ta:
            continue
        terms.extend(ta.get('hgvs', []))
        prot_hgvs = ta.get('proteinEffect', {}).get('hgvs')
        if prot_hgvs:
            terms.append(prot_hgvs)
    return '|'.join(terms) if terms else '-'


def _ucsc_url(chrom, pos, assembly):
    db = 'hg38' if assembly == 'GRCh38' else 'hg19'
    c = chrom if str(chrom).startswith('chr') else f'chr{chrom}'
    return f'https://genome.ucsc.edu/cgi-bin/hgTracks?db={db}&position={c}:{pos}-{pos}'


def _grch37_coords(data, vrs_digest):
    """Return a list of genomic_coordinates rows for GRCh37 genomicAlleles.

    Each row is a dict ready for INSERT into genomic_coordinates.
    The NC_ accession HGVS is preferred; falls back to the first hgvs entry.
    """
    rows = []
    for ga in data.get('genomicAlleles', []):
        if ga.get('referenceGenome') != 'GRCh37':
            continue
        coords = ga.get('coordinates', [{}])[0]
        chrom = ga.get('chromosome', '-')
        pos   = str(coords.get('end', '-'))
        ref   = coords.get('referenceAllele', '-')
        alt   = coords.get('allele', '-')
        hgvs_list = ga.get('hgvs', [])
        hgvs  = next((h for h in hgvs_list if h.startswith('NC_')), hgvs_list[0] if hgvs_list else '-')
        rows.append({
            'VRS_Digest_id': vrs_digest,
            'assembly':      'GRCh37',
            'hgvs':          hgvs,
            'genome_browser_url': _ucsc_url(chrom, pos, 'GRCh37'),
            'chr':  chrom,
            'pos':  pos,
            'ref':  ref,
            'alt':  alt,
        })
    return rows


def _extract(data, vrs_digest):
    """Return (variant_update_tuple, [coord_rows]) from a ClinGen response."""
    at_id = data.get('@id', '')
    ca_id = at_id.rsplit('/', 1)[-1] if at_id else None

    titles = data.get('communityStandardTitle')
    title = titles[0] if titles else None

    hgvs_cdna, ensembl_cdna, hgvs_protein, ensembl_protein = _mane_select(data)
    synonyms = _synonyms(data)
    coord_rows = _grch37_coords(data, vrs_digest)

    return ca_id, title, hgvs_cdna, ensembl_cdna, hgvs_protein, ensembl_protein, synonyms, coord_rows


@click.command()
@click.option('--db-url', default='postgresql://postgres:postgres@localhost/storage.pg',
              envvar='PIPELINE_DB_URL', show_default=True,
              help='PostgreSQL connection URL for the pipeline DB')
@click.option('--schema', default='pipeline', show_default=True,
              help='Schema containing variant and genomic_coordinates')
@click.option('--overwrite', is_flag=True, default=False,
              help='Update variants that already have a CA_ID')
def main(db_url, schema, overwrite):
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
                    JOIN variant v ON v."VRS_Digest" = gc."VRS_Digest_id"
                    WHERE gc.assembly = 'GRCh38'
                      AND gc.hgvs IS NOT NULL AND gc.hgvs <> ''
                      AND (v."CA_ID" IS NULL OR v."CA_ID" = '')
                """)
            rows = cur.fetchall()

        total = len(rows)
        print(f'Querying ClinGen for {total} variants ...')

        variant_updates = []
        coord_inserts   = []
        not_found = 0
        failed    = 0

        session = requests.Session()
        session.headers['Accept'] = 'application/json'

        for i, (vrs_digest, hgvs) in enumerate(rows, 1):
            if i % 500 == 0 or i == total:
                print(f'  {i}/{total}  ({len(variant_updates)} variant updates, '
                      f'{len(coord_inserts)} coord inserts, '
                      f'{not_found} not found, {failed} errors)')
            try:
                data = _query(session, hgvs)
            except Exception as e:
                logging.warning('Failed for %s (%s): %s', vrs_digest, hgvs, e)
                failed += 1
                continue

            if data is None:
                not_found += 1
                continue

            ca_id, title, hgvs_cdna, ensembl_cdna, hgvs_protein, ensembl_protein, synonyms, coord_rows = \
                _extract(data, vrs_digest)

            if ca_id:
                ref_seq = hgvs_cdna.split(':')[0] if hgvs_cdna and ':' in hgvs_cdna else '-'
                variant_updates.append((
                    ca_id,
                    title or '-',
                    hgvs_cdna or '-',
                    ensembl_cdna or '-',
                    hgvs_protein or '-',
                    ensembl_protein or '-',
                    synonyms,
                    ref_seq,
                    vrs_digest,
                ))

            coord_inserts.extend(coord_rows)

            time.sleep(REQUEST_DELAY)

        print(f'Writing {len(variant_updates)} variant updates and '
              f'{len(coord_inserts)} GRCh37 coordinate inserts ...')

        with conn.cursor() as cur:
            for i in range(0, len(variant_updates), BATCH_SIZE):
                batch = variant_updates[i : i + BATCH_SIZE]
                cur.executemany(
                    """UPDATE variant
                       SET "CA_ID" = %s, title = %s, "HGVS_cDNA" = %s, ensembl_cdna = %s,
                           "HGVS_Protein" = %s, ensembl_protein = %s, "Synonyms" = %s,
                           "Reference_Sequence" = %s
                       WHERE "VRS_Digest" = %s""",
                    batch,
                )

            for i in range(0, len(coord_inserts), BATCH_SIZE):
                batch = coord_inserts[i : i + BATCH_SIZE]
                psycopg2.extras.execute_values(
                    cur,
                    """INSERT INTO genomic_coordinates
                           ("VRS_Digest_id", assembly, hgvs, genome_browser_url, chr, pos, ref, alt)
                       VALUES %s
                       ON CONFLICT ("VRS_Digest_id", assembly) DO NOTHING""",
                    [(r['VRS_Digest_id'], r['assembly'], r['hgvs'],
                      r['genome_browser_url'], r['chr'], r['pos'], r['ref'], r['alt'])
                     for r in batch],
                )

        conn.commit()
        print(f'Done. Variants updated: {len(variant_updates)}, '
              f'GRCh37 coords inserted: {len(coord_inserts)}, '
              f'not found: {not_found}, errors: {failed}.')
    finally:
        conn.close()


if __name__ == '__main__':
    main()
