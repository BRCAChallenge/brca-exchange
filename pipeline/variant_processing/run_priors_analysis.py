"""
Populate analysis_priors from calcVarPriors output.

Processes substitution variants only — the splicing prior calculator
handles only SNVs. Uses varType from analysis_vep (populated by
run_vep_analysis.py) to filter variants without re-querying the VEP
server. Writes a temp input TSV, runs calcVarPriors.py calc as a
subprocess (from the splicing directory), parses the output, and upserts
the UI-visible priors columns into analysis_priors.
"""

import csv
import os
import subprocess
import tempfile

import click
import psycopg2
import psycopg2.extras

DB_BATCH = 500

_SPLICING_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    'splicing',
)

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
@click.option('--genome', default='/references/hg38.fa', show_default=True,
              help='Path to hg38 reference FASTA')
@click.option('--processes', default=8, show_default=True,
              help='Number of parallel calcVarPriors workers')
@click.option('--overwrite', is_flag=True, default=False,
              help='Re-score variants already in analysis_priors')
def main(db_url, schema, genome, processes, overwrite):
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
                      AND ap."VRS_Digest" IS NULL
                """)
            db_rows = cur.fetchall()
    finally:
        conn.close()

    snvs = []
    for vrs, gene, refseq, hgvs_cdna, chr_, pos, ref, alt in db_rows:
        chr_bare = chr_[3:] if chr_.startswith('chr') else chr_
        snvs.append((vrs, gene, refseq, hgvs_cdna, chr_bare, pos, ref, alt))

    print(f'Processing {len(snvs)} SNV variants ...')
    if not snvs:
        print('Nothing to do.')
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        in_tsv = os.path.join(tmpdir, 'variants_in.tsv')
        out_tsv = os.path.join(tmpdir, 'variants_out.tsv')

        with open(in_tsv, 'w', newline='') as fh:
            writer = csv.writer(fh, delimiter='\t')
            writer.writerow(['Chr', 'Pos', 'Ref', 'Alt', 'Reference_Sequence',
                             'HGVS_cDNA', 'Gene_Symbol', 'Hg38_Start', 'Hg38_End', 'VRS_Digest'])
            for vrs, gene, refseq, hgvs_cdna, chr_bare, pos, ref, alt in snvs:
                writer.writerow([chr_bare, pos, ref, alt, refseq,
                                  hgvs_cdna, gene, pos, pos, vrs])

        cmd = [
            'python', 'calcVarPriors.py',
            '--genome', genome,
            '--processes', str(processes),
            'calc', in_tsv, out_tsv,
        ]
        print(f'Running calcVarPriors.py ...')
        result = subprocess.run(cmd, cwd=_SPLICING_DIR)
        if result.returncode != 0:
            raise SystemExit(f'calcVarPriors.py exited with code {result.returncode}')

        print(f'Parsing output ...')
        inserts = []
        with open(out_tsv, newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                if row.get('varType') != 'substitution':
                    continue
                vrs = row['VRS_Digest']
                values = tuple(_blank_to_none(row.get(col)) for col in _PRIORS_COLS)
                inserts.append((vrs,) + values)

    print(f'Writing {len(inserts)} rows to analysis_priors ...')
    col_list = ', '.join(f'"{c}"' for c in _PRIORS_COLS)
    update_set = ', '.join(f'"{c}" = EXCLUDED."{c}"' for c in _PRIORS_COLS)

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
