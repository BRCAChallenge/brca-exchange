"""
Populate analysis_priors from calcVarPriors.

Processes substitution variants only — the splicing prior calculator
handles only SNVs. Uses variant_type from analysis_vep (populated by
run_vep_analysis.py) to filter variants without re-querying the VEP
server. Imports calc_one from calcVarPriors directly and upserts the
UI-visible priors columns into analysis_priors.

Variants whose Reference_Sequence:HGVS_cDNA (version-normalized) appears
in splicingfilter/blacklisted_vars.txt are skipped; their priors columns
are left NULL in analysis_priors.
"""

import io
import multiprocessing
import os
import sys

import click
import psycopg2
import psycopg2.extras
import pyhgvs.utils as pyhgvs_utils

_DIR = os.path.dirname(os.path.abspath(__file__))
_PRIORS_DIR = os.path.join(_DIR, 'priors')
sys.path.insert(0, _PRIORS_DIR)

import calcVarPriors
from calc_priors.constants import BRCA1_RefSeq, BRCA2_RefSeq

DB_BATCH = 500
_BLACKLIST_PATH = os.path.join(_DIR, '..', 'splicingfilter', 'blacklisted_vars.txt')


def _load_blacklist(path):
    """Return a set of version-stripped 'ACCESSION:cDNA' keys, e.g. 'NM_000059:c.551T>C'."""
    blacklist = set()
    try:
        with open(path) as f:
            for line in f:
                entry = line.strip()
                if not entry:
                    continue
                if ':' in entry:
                    acc, cdna = entry.split(':', 1)
                    # strip version suffix so NM_000059.3 matches NM_000059.4
                    acc_base = acc.rsplit('.', 1)[0]
                    blacklist.add(f'{acc_base}:{cdna}')
    except FileNotFoundError:
        print(f'Warning: blacklist file not found at {path}; no variants will be skipped.')
    return blacklist


def _blacklist_key(refseq, hgvs_cdna):
    """Return version-stripped lookup key for a variant."""
    acc_base = refseq.rsplit('.', 1)[0] if refseq else refseq
    return f'{acc_base}:{hgvs_cdna}'

# UI-visible priors columns stored in analysis_priors
_PRIORS_COLS = [
    'varLoc',
    'proteinPrior',
    'refDonorPrior',
    'refRefDonorMES', 'refRefDonorZ', 'altRefDonorMES', 'altRefDonorZ',
    'refRefDonorSeq', 'altRefDonorSeq',
    'deNovoDonorPrior',
    'refDeNovoDonorMES', 'refDeNovoDonorZ', 'altDeNovoDonorMES', 'altDeNovoDonorZ',
    'refDeNovoDonorSeq', 'altDeNovoDonorSeq',
    'deNovoDonorGenomicSplicePos', 'deNovoDonorTranscriptSplicePos',
    'closestDonorGenomicSplicePos', 'closestDonorTranscriptSplicePos',
    'closestDonorRefMES', 'closestDonorRefZ', 'closestDonorRefSeq',
    'closestDonorAltMES', 'closestDonorAltZ', 'closestDonorAltSeq',
    'refAccPrior',
    'refRefAccMES', 'refRefAccZ', 'altRefAccMES', 'altRefAccZ',
    'refRefAccSeq', 'altRefAccSeq',
]


def _blank_to_none(val):
    return None if val in ('-', '', None) else val


@click.command()
@click.option('--db-url', default='postgresql://postgres:postgres@localhost/storage.pg',
              envvar='PIPELINE_DB_URL', show_default=True)
@click.option('--schema', default='pipeline', show_default=True)
@click.option('--processes', default=8, show_default=True,
              help='Number of parallel calcVarPriors workers')
@click.option('--overwrite', is_flag=True, default=False,
              help='Re-score variants already in analysis_priors')
@click.option('--limit', default=0, show_default=True,
              help='Process only this many variants (0 = all)')
def main(db_url, schema, processes, overwrite, limit):
    transcripts_path = os.path.join(_PRIORS_DIR, 'refseq_annotation.hg38.gp')
    with open(transcripts_path) as f:
        # read_transcripts requires 16-column genePredExt format; skip shorter lines
        filtered = io.StringIO(''.join(
            line for line in f if line.startswith('#') or len(line.split('\t')) == 16
        ))
    transcripts = pyhgvs_utils.read_transcripts(filtered)
    calcVarPriors.brca1Transcript = transcripts.get(BRCA1_RefSeq)
    calcVarPriors.brca2Transcript = transcripts.get(BRCA2_RefSeq)

    conn = psycopg2.connect(db_url, options=f'-c search_path={schema}')
    try:
        with conn.cursor() as cur:
            if overwrite:
                cur.execute("""
                    SELECT v."VRS_Digest", v."Gene_Symbol", v."Reference_Sequence", v."HGVS_cDNA",
                           gc.chr, gc.pos, gc.ref, gc.alt
                    FROM variant v
                    JOIN genomic_coordinates gc ON gc."VRS_Digest_id" = v."VRS_Digest"
                    JOIN analysis_vep av ON av."VRS_Digest" = v."VRS_Digest"
                    WHERE gc.assembly = 'GRCh38'
                      AND av.variant_type = 'substitution'
                      AND v."Gene_Symbol" IN ('BRCA1', 'BRCA2')
                """)
            else:
                cur.execute("""
                    SELECT v."VRS_Digest", v."Gene_Symbol", v."Reference_Sequence", v."HGVS_cDNA",
                           gc.chr, gc.pos, gc.ref, gc.alt
                    FROM variant v
                    JOIN genomic_coordinates gc ON gc."VRS_Digest_id" = v."VRS_Digest"
                    JOIN analysis_vep av ON av."VRS_Digest" = v."VRS_Digest"
                    LEFT JOIN analysis_priors ap ON ap."VRS_Digest" = v."VRS_Digest"
                    WHERE gc.assembly = 'GRCh38'
                      AND av.variant_type = 'substitution'
                      AND v."Gene_Symbol" IN ('BRCA1', 'BRCA2')
                      AND ap."VRS_Digest" IS NULL
                """)
            db_rows = cur.fetchall()
    finally:
        conn.close()

    blacklist = _load_blacklist(_BLACKLIST_PATH)

    variants = []
    blacklisted_count = 0
    for vrs, gene, refseq, hgvs_cdna, chr_, pos, ref, alt in db_rows:
        if blacklist and _blacklist_key(refseq, hgvs_cdna) in blacklist:
            blacklisted_count += 1
            continue
        chr_bare = chr_[3:] if chr_.startswith('chr') else chr_
        variants.append({
            'VRS_Digest': vrs,
            'Gene_Symbol': gene,
            'Reference_Sequence': refseq,
            'HGVS_cDNA': hgvs_cdna,
            'Chr': chr_bare,
            'Pos': pos,
            'Ref': ref,
            'Alt': alt,
            'Hg38_Start': pos,
            'Hg38_End': pos,
        })

    if blacklisted_count:
        print(f'Skipped {blacklisted_count} blacklisted variants.')
    if limit > 0:
        variants = variants[:limit]
    print(f'Processing {len(variants)} substitution variants ...')
    if not variants:
        print('Nothing to do.')
        return

    if processes > 1:
        pool = multiprocessing.Pool(processes)
        try:
            rows = pool.map(calcVarPriors.calc_one, variants)
        finally:
            pool.close()
            pool.join()
    else:
        rows = [calcVarPriors.calc_one(v) for v in variants]

    inserts = []
    for row in rows:
        vrs = row['VRS_Digest']
        values = tuple(_blank_to_none(row.get(col)) for col in _PRIORS_COLS)
        inserts.append((vrs,) + values)

    col_list = ', '.join(f'"{c}"' for c in _PRIORS_COLS)
    update_set = ', '.join(f'"{c}" = EXCLUDED."{c}"' for c in _PRIORS_COLS)

    print(f'Writing {len(inserts)} rows to analysis_priors ...')
    conn = psycopg2.connect(db_url, options=f'-c search_path={schema}')
    try:
        with conn.cursor() as cur:
            for i in range(0, len(inserts), DB_BATCH):
                psycopg2.extras.execute_values(
                    cur,
                    f"""INSERT INTO analysis_priors ("VRS_Digest", {col_list})
                        VALUES %s
                        ON CONFLICT ("VRS_Digest") DO UPDATE SET {update_set}""",
                    inserts[i: i + DB_BATCH],
                )
        conn.commit()
        print(f'Done. Inserted/updated {len(inserts)} rows.')
    finally:
        conn.close()


if __name__ == '__main__':
    main()
