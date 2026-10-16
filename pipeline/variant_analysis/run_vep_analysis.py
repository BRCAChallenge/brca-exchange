"""
Run VEP on all variants and populate analysis_vep.

For each variant with a GRCh38 genomic HGVS in variant_genomic_coordinates,
queries the local VEP server in batches (against its RefSeq-cache endpoint)
and stores variant_class, varType, consequences, hgvsp, and (for variants
that may introduce a premature stop codon) ptc_genomic_pos in the
analysis_vep table.

varType is the structural classification of the variant (substitution,
insertion, deletion, delins, other) derived from ref/alt — the same
classification used by the splicing priors calculator.

consequences is the comma-separated list of VEP consequence_terms
(e.g. "missense_variant"), and hgvsp is the VEP HGVSp notation (e.g.
"NP_009225.1:p.Gln356GlufsTer9"), both taken from the transcript_consequences
entry matching the variant's Reference_Sequence RefSeq transcript accession
(variant.Reference_Sequence, e.g. NM_007294.4) — not necessarily the gene's
canonical transcript, since Reference_Sequence varies per variant.

ptc_genomic_pos is the GRCh38 1-based genomic position of the first base of
the premature stop codon, for variants whose matched consequence_terms
intersect PTC_CONSEQUENCE_TERMS ({'stop_gained', 'frameshift_variant'}).
VEP's own core annotation (--hgvs) already computes the stop's position in
protein coordinates as part of hgvsp — for a frameshift it translates the
shifted reading frame codon-by-codon until an in-frame stop and reports the
offset as "TerN" (e.g. p.Gln356GlufsTer9 = new stop 9 residues downstream of
residue 356, i.e. protein position 364). This is parsed from hgvsp rather
than re-derived. What VEP does not do is project that protein position back
to a genomic coordinate through the transcript's exon structure — that step
uses the biocommons hgvs library against a UTA database instance, the same
approach used elsewhere in this pipeline (see run_priors_analysis.py). Set
UTA_DB_URL to point at a local instance for best performance, e.g.:
  UTA_DB_URL=postgresql://anonymous@localhost:50828/uta/uta_20241220

Variants that already have a row in analysis_vep are skipped unless
--overwrite is set.
"""

import logging
import os
import re

import click
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.exceptions
import hgvs.parser
import psycopg2
import psycopg2.extras
import requests

BATCH_SIZE   = 200   # variants per VEP call
DB_BATCH     = 500   # rows per INSERT batch

log = logging.getLogger(__name__)

_ACCEPTABLE = set('ACGTNRY')

# Mirrors PTC_CONSEQUENCE_TERMS in django/data/models.py — kept in sync by hand
# since this script runs standalone (psycopg2, no Django app registry) and
# can't import the Django model directly.
PTC_CONSEQUENCE_TERMS = {'stop_gained', 'frameshift_variant'}

# VEP hgvsp protein change, e.g. "Gln356GlufsTer9" (frameshift) or "Trp964Ter"
# (direct nonsense). Ter offset may be '?' when VEP can't find an in-frame
# stop before the transcript ends.
_FS_RE   = re.compile(r'p\.\(?[A-Za-z]{3}(\d+)[A-Za-z]{3}fs(?:Ter|\*)(\d+|\?)\)?$')
_STOP_RE = re.compile(r'p\.\(?[A-Za-z]{3}(\d+)(?:Ter|\*)\)?$')

_hp = hgvs.parser.Parser()
_mapper = None  # initialised lazily so UTA_DB_URL can be set at runtime


def _get_mapper():
    global _mapper
    if _mapper is None:
        _mapper = hgvs.assemblymapper.AssemblyMapper(
            hgvs.dataproviders.uta.connect(), assembly_name='GRCh38',
            alt_aln_method='splign', prevalidation_level=None)
    return _mapper


def _ptc_stop_protein_pos(hgvsp):
    """Given a VEP hgvsp string (e.g. 'NP_009225.1:p.Gln356GlufsTer9' or
    'NP_009225.1:p.Trp964Ter'), return the protein position (1-based) of the
    first residue of the resulting premature stop codon, or None if hgvsp is
    missing, doesn't describe a stop, or the frameshift's stop position is
    unresolved ('Ter?').
    """
    if not hgvsp or ':' not in hgvsp:
        return None
    p_part = hgvsp.split(':', 1)[1]
    m = _FS_RE.match(p_part)
    if m:
        offset = m.group(2)
        if offset == '?':
            return None
        return int(m.group(1)) + int(offset) - 1
    m = _STOP_RE.match(p_part)
    if m:
        return int(m.group(1))
    return None


def _cds_pos_to_genomic_pos(transcript, cds_pos):
    """Map a CDS position (1-based, c. numbering from ATG=1) on transcript to
    a GRCh38 genomic position, via a dummy substitution — prevalidation is
    disabled on the mapper, so the placeholder ref/alt bases are never
    checked against real sequence, only the position is used. Mirrors the
    same trick used in calc_priors/verify.py's convertGenomicPosToTranscriptPos,
    in the reverse direction.
    """
    c_var = _hp.parse_hgvs_variant(f'{transcript}:c.{cds_pos}A>T')
    g_var = _get_mapper().c_to_g(c_var)
    return int(g_var.posedit.pos.start.base)


def _get_var_type(ref, alt):
    """Classify a variant by structural type, matching calcVarPriors getVarType logic."""
    ref = ref.upper() if ref else ''
    alt = alt.upper() if alt else ''
    ok = (bool(ref) and all(c in _ACCEPTABLE for c in ref)
          and bool(alt) and all(c in _ACCEPTABLE for c in alt))
    if not ok:
        return 'other'
    if len(ref) == len(alt):
        return 'substitution' if len(ref) == 1 else 'delins'
    if len(ref) > len(alt):
        return 'deletion' if len(alt) == 1 and alt == ref[0] else 'delins'
    # len(ref) < len(alt)
    return 'insertion' if len(ref) == 1 and ref == alt[0] else 'delins'


def _query_vep_batch(session, vep_url, hgvs_list):
    """POST a batch of HGVS strings to the VEP server's RefSeq-cache endpoint.
    Returns dict {hgvs: [records]} or {hgvs: {'error': ...}}.
    """
    resp = session.post(f'{vep_url}/vep/hgvs/batch/refseq',
                        json=hgvs_list, timeout=600)
    resp.raise_for_status()
    return resp.json()


def _match_transcript_consequence(records, ref_seq):
    """Find the transcript_consequences entry for ref_seq (e.g. 'NM_007294.4')
    among VEP result records, and return that entry (dict), or None if no
    matching transcript is found.

    Falls back to a version-insensitive match (base accession only) since the
    RefSeq version in the VEP cache can lag the version stored on the variant.
    """
    if not ref_seq or ref_seq == '-' or not records:
        return None
    base = ref_seq.split('.')[0]
    fallback = None
    for record in records:
        for tc in record.get('transcript_consequences', []):
            tid = tc.get('transcript_id', '')
            if tid == ref_seq:
                return tc
            if fallback is None and tid.split('.')[0] == base:
                fallback = tc
    return fallback


@click.command()
@click.option('--db-url', default='postgresql://postgres:postgres@localhost/storage.pg',
              envvar='PIPELINE_DB_URL', show_default=True)
@click.option('--schema', default='pipeline', show_default=True)
@click.option('--vep-url', default='http://localhost:8888', show_default=True,
              envvar='VEP_SERVER_URL', help='Base URL of the local VEP REST server')
@click.option('--overwrite', is_flag=True, default=False,
              help='Re-annotate variants that already have a row in analysis_vep')
@click.option('--uta-db-url', default=None, envvar='UTA_DB_URL', show_default=True,
              help='UTA database URL for transcript coordinate mapping (PTC genomic '
                   'position). Defaults to UTA_DB_URL env var or the biocommons public '
                   'instance. Use a local instance for best performance, e.g. '
                   'postgresql://anonymous@localhost:50828/uta/uta_20241220')
def main(db_url, schema, vep_url, overwrite, uta_db_url):
    logging.basicConfig(level=logging.WARNING, format='%(levelname)s %(message)s')

    if uta_db_url:
        os.environ['UTA_DB_URL'] = uta_db_url

    conn = psycopg2.connect(db_url, options=f'-c search_path={schema}')
    try:
        with conn.cursor() as cur:
            if overwrite:
                cur.execute("""
                    SELECT gc."VRS_Digest", gc.hgvs, gc.ref, gc.alt, v."Reference_Sequence"
                    FROM variant_genomic_coordinates gc
                    JOIN variant v ON v."VRS_Digest" = gc."VRS_Digest"
                    WHERE gc.assembly = 'GRCh38'
                      AND gc.hgvs IS NOT NULL AND gc.hgvs <> ''
                """)
            else:
                cur.execute("""
                    SELECT gc."VRS_Digest", gc.hgvs, gc.ref, gc.alt, v."Reference_Sequence"
                    FROM variant_genomic_coordinates gc
                    JOIN variant v ON v."VRS_Digest" = gc."VRS_Digest"
                    LEFT JOIN analysis_vep av ON av."VRS_Digest" = gc."VRS_Digest"
                    WHERE gc.assembly = 'GRCh38'
                      AND gc.hgvs IS NOT NULL AND gc.hgvs <> ''
                      AND av."VRS_Digest" IS NULL
                """)
            rows = cur.fetchall()

        # Build lookups so varType/consequences are available when VEP results come back
        ref_alt = {r[0]: (r[2], r[3]) for r in rows}   # VRS_Digest -> (ref, alt)
        ref_seq = {r[0]: r[4] for r in rows}           # VRS_Digest -> Reference_Sequence

        total = len(rows)
        print(f'Running VEP on {total} variants in batches of {BATCH_SIZE} ...')

        session = requests.Session()
        inserts = []
        errors  = 0

        for batch_start in range(0, total, BATCH_SIZE):
            batch = rows[batch_start : batch_start + BATCH_SIZE]
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

            for hgvs_str, records in results.items():
                vrs = hgvs_map.get(hgvs_str)
                if not vrs:
                    continue
                if isinstance(records, dict) and 'error' in records:
                    log.warning('VEP error for %s: %s', hgvs_str, records['error'])
                    errors += 1
                    continue
                variant_class = records[0].get('variant_class') if records else None
                ref, alt = ref_alt.get(vrs, (None, None))
                var_type = _get_var_type(ref, alt) if ref is not None else None
                tc = _match_transcript_consequence(records, ref_seq.get(vrs))
                consequences = ','.join(tc.get('consequence_terms', [])) if tc else None
                hgvsp = tc.get('hgvsp') if tc else None

                ptc_genomic_pos = None
                if tc and set(tc.get('consequence_terms', [])) & PTC_CONSEQUENCE_TERMS:
                    stop_pos = _ptc_stop_protein_pos(hgvsp)
                    if stop_pos:
                        try:
                            ptc_genomic_pos = _cds_pos_to_genomic_pos(
                                ref_seq[vrs], 3 * stop_pos - 2)
                        except hgvs.exceptions.HGVSError as e:
                            log.warning('PTC genomic mapping failed for %s: %s', hgvs_str, e)

                inserts.append((vrs, variant_class, var_type, consequences, hgvsp, ptc_genomic_pos))

        print(f'Writing {len(inserts)} rows to analysis_vep ...')
        with conn.cursor() as cur:
            for i in range(0, len(inserts), DB_BATCH):
                psycopg2.extras.execute_values(
                    cur,
                    """INSERT INTO analysis_vep
                           ("VRS_Digest", variant_class, variant_type, consequences,
                            hgvsp, ptc_genomic_pos)
                       VALUES %s
                       ON CONFLICT ("VRS_Digest") DO UPDATE
                         SET variant_class    = EXCLUDED.variant_class,
                             variant_type     = EXCLUDED.variant_type,
                             consequences     = EXCLUDED.consequences,
                             hgvsp            = EXCLUDED.hgvsp,
                             ptc_genomic_pos  = EXCLUDED.ptc_genomic_pos""",
                    inserts[i : i + DB_BATCH],
                )
        conn.commit()
        print(f'Done. Inserted/updated {len(inserts)}, errors {errors}.')
    finally:
        conn.close()


if __name__ == '__main__':
    main()
