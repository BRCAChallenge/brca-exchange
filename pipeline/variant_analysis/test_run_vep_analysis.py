"""
Unit tests for run_vep_analysis.py.

Pure-logic tests only -- no database or VEP/UTA network connection is needed,
since _ptc_stop_protein_pos and _match_transcript_consequence are plain
functions. main() itself (the click entrypoint that talks to PostgreSQL and
VEP) and _cds_pos_to_genomic_pos (talks to UTA) are intentionally not
exercised here.

Run with:
    pytest variant_analysis/test_run_vep_analysis.py
"""

import pytest

import run_vep_analysis as vep


# ---------------------------------------------------------------------------
# _ptc_stop_protein_pos
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("hgvsp,expected", [
    # frameshift: new stop is 9 residues downstream of residue 356 -> 364
    ("NP_009225.1:p.Gln356GlufsTer9", 364),
    ("NP_009225.1:p.Gln356Glufs*9", 364),
    ("NP_009225.1:p.(Gln356GlufsTer9)", 364),
    # direct nonsense: stop replaces residue 964 itself
    ("NP_009225.1:p.Trp964Ter", 964),
    ("NP_009225.1:p.Trp964*", 964),
    # unresolved frameshift length -> can't compute a position
    ("NP_009225.1:p.Ile327ArgfsTer?", None),
    # not a stop at all (missense) -> no PTC
    ("NP_009225.1:p.Gly12Val", None),
    # missing/malformed hgvsp
    (None, None),
    ("", None),
    ("garbage", None),
])
def test_ptc_stop_protein_pos(hgvsp, expected):
    assert vep._ptc_stop_protein_pos(hgvsp) == expected


# ---------------------------------------------------------------------------
# _match_transcript_consequence
# ---------------------------------------------------------------------------

def _records(*tcs):
    return [{'transcript_consequences': list(tcs)}]


def test_match_transcript_consequence_exact_match():
    tc_match = {'transcript_id': 'NM_007294.4', 'consequence_terms': ['missense_variant']}
    tc_other = {'transcript_id': 'NM_000059.4', 'consequence_terms': ['missense_variant']}
    records = _records(tc_other, tc_match)
    assert vep._match_transcript_consequence(records, 'NM_007294.4') is tc_match


def test_match_transcript_consequence_version_insensitive_fallback():
    tc_older_version = {'transcript_id': 'NM_007294.3', 'consequence_terms': ['stop_gained']}
    records = _records(tc_older_version)
    assert vep._match_transcript_consequence(records, 'NM_007294.4') is tc_older_version


def test_match_transcript_consequence_no_match():
    tc_other = {'transcript_id': 'NM_000059.4', 'consequence_terms': ['missense_variant']}
    records = _records(tc_other)
    assert vep._match_transcript_consequence(records, 'NM_007294.4') is None


@pytest.mark.parametrize("ref_seq,records", [
    (None, _records({'transcript_id': 'NM_007294.4'})),
    ('-', _records({'transcript_id': 'NM_007294.4'})),
    ('NM_007294.4', []),
])
def test_match_transcript_consequence_missing_inputs(ref_seq, records):
    assert vep._match_transcript_consequence(records, ref_seq) is None
