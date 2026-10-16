"""
Tests for pseudonym_generator HGVS_cDNA fill-in behaviour.

Focus: ensure_mane_transcript_cdna and the surrounding process_rows flow,
specifically the cases where HGVS_cDNA is not supplied by the source VCF
and must be derived from SeqRepo / ClinGen.
"""

import copy
import os
import re
from unittest.mock import MagicMock, patch, call
import pytest
import hgvs.exceptions

from variant_assembly.pseudonym_generator import (
    ensure_mane_transcript_cdna,
    _init_hgvs_tools,
    HGVS_CDNA_COL,
    PYHGVS_CDNA_COL,
    PYHGVS_GENOMIC_COORDINATE_38_COL,
    REFERENCE_SEQUENCE_COL,
)

MANE = 'NM_007294.4'


def _row(hgvs_cdna='-', pyhgvs_cdna=None, ref_seq='-',
         genomic_38='NC_000017.11:g.43094692G>A'):
    """Return a minimal row dict as produced by _load_rows_from_db / process_rows."""
    return {
        HGVS_CDNA_COL:                   hgvs_cdna,
        PYHGVS_CDNA_COL:                  pyhgvs_cdna,
        REFERENCE_SEQUENCE_COL:           ref_seq,
        PYHGVS_GENOMIC_COORDINATE_38_COL: genomic_38,
    }


def _mock_hgvs_proc(parsed_obj=None):
    proc = MagicMock()
    proc.hgvs_parser.parse.return_value = parsed_obj or MagicMock()
    return proc


# ---------------------------------------------------------------------------
# ensure_mane_transcript_cdna
# ---------------------------------------------------------------------------

class TestEnsureManeTranscriptCdna:

    def test_seqrepo_fills_missing_hgvs(self):
        """When HGVS_cDNA is '-' and SeqRepo succeeds, the row gets the cDNA value."""
        row = _row(hgvs_cdna='-', pyhgvs_cdna=None)
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()

        cdna_result = MagicMock()
        cdna_result.__str__ = lambda self: f'{MANE}:c.5266dupC'
        am38.g_to_c.return_value = cdna_result
        normalizer.normalize.return_value = cdna_result

        ensure_mane_transcript_cdna(row, MANE, hgvs_proc, normalizer, am38)

        assert row[HGVS_CDNA_COL] == 'c.5266dupC'
        assert row[REFERENCE_SEQUENCE_COL] == MANE

    def test_seqrepo_fills_missing_hgvs_brca2(self):
        """Same check for a BRCA2 MANE transcript."""
        mane2 = 'NM_000059.4'
        row = _row(hgvs_cdna='-', pyhgvs_cdna=None,
                   genomic_38='NC_000013.11:g.32340300A>G')
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()

        cdna_result = MagicMock()
        cdna_result.__str__ = lambda self: f'{mane2}:c.1114A>G'
        am38.g_to_c.return_value = cdna_result
        normalizer.normalize.return_value = cdna_result

        ensure_mane_transcript_cdna(row, mane2, hgvs_proc, normalizer, am38)

        assert row[HGVS_CDNA_COL] == 'c.1114A>G'

    def test_clingen_cdna_accepted_without_seqrepo_call(self):
        """When ClinGen already returned the MANE cDNA, SeqRepo is not called."""
        pyhgvs_cdna = f'{MANE}:c.5266dupC'
        row = _row(hgvs_cdna='-', pyhgvs_cdna=pyhgvs_cdna, ref_seq=MANE)
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()

        ensure_mane_transcript_cdna(row, MANE, hgvs_proc, normalizer, am38)

        # SeqRepo must not have been called; value comes from ClinGen result
        am38.g_to_c.assert_not_called()
        assert row[HGVS_CDNA_COL] == 'c.5266dupC'

    def test_seqrepo_failure_preserves_existing_source_hgvs(self):
        """When SeqRepo fails but the row already has a VCF-sourced HGVS_cDNA,
        the existing value is kept (the key bug fix)."""
        existing = 'c.-40+1G>A'
        row = _row(hgvs_cdna=existing, pyhgvs_cdna=None)
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()
        am38.g_to_c.side_effect = hgvs.exceptions.HGVSInvalidIntervalError('out of range')

        ensure_mane_transcript_cdna(row, MANE, hgvs_proc, normalizer, am38)

        assert row[HGVS_CDNA_COL] == existing, (
            "SeqRepo failure must not overwrite a pre-existing HGVS_cDNA value"
        )

    def test_seqrepo_failure_no_existing_hgvs_stays_dash(self):
        """When SeqRepo fails and there was no pre-existing HGVS_cDNA,
        the value remains '-' (correct for gnomAD-only deep intronic variants)."""
        row = _row(hgvs_cdna='-', pyhgvs_cdna=None)
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()
        am38.g_to_c.side_effect = hgvs.exceptions.HGVSInvalidIntervalError('out of range')

        ensure_mane_transcript_cdna(row, MANE, hgvs_proc, normalizer, am38)

        assert row[HGVS_CDNA_COL] == '-'

    def test_seqrepo_unsupported_operation_preserves_existing(self):
        """HGVSUnsupportedOperationError is also caught and handled correctly."""
        existing = 'c.67+58A>C'
        row = _row(hgvs_cdna=existing, pyhgvs_cdna=None)
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()
        am38.g_to_c.side_effect = hgvs.exceptions.HGVSUnsupportedOperationError('unsupported')

        ensure_mane_transcript_cdna(row, MANE, hgvs_proc, normalizer, am38)

        assert row[HGVS_CDNA_COL] == existing

    def test_seqrepo_data_not_available_preserves_existing(self):
        """HGVSDataNotAvailableError is also caught and handled correctly."""
        existing = 'c.67+82C>G'
        row = _row(hgvs_cdna=existing, pyhgvs_cdna=None)
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()
        am38.g_to_c.side_effect = hgvs.exceptions.HGVSDataNotAvailableError('no data')

        ensure_mane_transcript_cdna(row, MANE, hgvs_proc, normalizer, am38)

        assert row[HGVS_CDNA_COL] == existing

    def test_seqrepo_overrides_non_mane_pyhgvs_cdna(self):
        """When ClinGen returns a cDNA on a non-MANE transcript, SeqRepo remaps it."""
        non_mane = 'NM_007294.3:c.5266dupC'
        row = _row(hgvs_cdna='-', pyhgvs_cdna=non_mane)
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()

        cdna_result = MagicMock()
        cdna_result.__str__ = lambda self: f'{MANE}:c.5266dupC'
        normalized_g = MagicMock()
        normalizer.normalize.return_value = normalized_g
        am38.g_to_c.return_value = cdna_result

        ensure_mane_transcript_cdna(row, MANE, hgvs_proc, normalizer, am38)

        am38.g_to_c.assert_called_once()
        assert row[HGVS_CDNA_COL] == 'c.5266dupC'
        assert row[REFERENCE_SEQUENCE_COL] == MANE

    def test_intronic_variant_gets_cdna_from_g_to_c(self):
        """Intronic variants: genomic is normalized first, g_to_c returns intronic cDNA
        notation (c.X+N or c.X-N), and cDNA normalization is NOT attempted."""
        row = _row(hgvs_cdna='-', pyhgvs_cdna=None,
                   genomic_38='NC_000017.11:g.43045803C>A')
        hgvs_proc = _mock_hgvs_proc()
        normalizer = MagicMock()
        am38 = MagicMock()

        normalized_g = MagicMock()
        normalizer.normalize.return_value = normalized_g

        intronic_cdna = MagicMock()
        intronic_cdna.__str__ = lambda self: f'{MANE}:c.5407-45C>A'
        am38.g_to_c.return_value = intronic_cdna

        ensure_mane_transcript_cdna(row, MANE, hgvs_proc, normalizer, am38)

        # Genomic normalization called once; cDNA normalization NOT called on the result
        normalizer.normalize.assert_called_once()
        am38.g_to_c.assert_called_once_with(normalized_g, MANE)
        assert row[HGVS_CDNA_COL] == 'c.5407-45C>A'
        assert row[REFERENCE_SEQUENCE_COL] == MANE


# ---------------------------------------------------------------------------
# Integration tests — require real SeqRepo/UTA
# Run with: HGVS_SEQREPO_DIR=/data/seqrepo/latest
#           UTA_DB_URL=postgresql://anonymous@localhost:50828/uta/uta_20241220
#           pytest -m integration
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def hgvs_tools():
    """Initialize real SeqRepo/UTA HGVS tools once per module."""
    os.environ.setdefault('HGVS_SEQREPO_DIR', '/data/seqrepo/latest')
    os.environ.setdefault('UTA_DB_URL',
                          'postgresql://anonymous@localhost:50828/uta/uta_20241220')
    hp, normalizer, am38, _am37, hgvs_proc, _lo = _init_hgvs_tools()
    return hgvs_proc, normalizer, am38


@pytest.mark.integration
class TestEnsureManeTranscriptCdnaIntegration:
    """Integration tests using real SeqRepo/UTA.

    Variant coordinates come from gnomAD-only rows in the pipeline schema that
    had HGVS_cDNA = '-' before the normalize-genomic-first fix.
    """

    BRCA1_MANE = 'NM_007294.4'
    BRCA2_MANE = 'NM_000059.4'

    def _run(self, genomic_38, mane, hgvs_tools):
        hgvs_proc, normalizer, am38 = hgvs_tools
        row = {
            HGVS_CDNA_COL:                   '-',
            PYHGVS_CDNA_COL:                  None,
            REFERENCE_SEQUENCE_COL:           '-',
            PYHGVS_GENOMIC_COORDINATE_38_COL: genomic_38,
        }
        ensure_mane_transcript_cdna(row, mane, hgvs_proc, normalizer, am38)
        return row

    def test_brca1_intronic_snv(self, hgvs_tools):
        """gnomAD-only intronic BRCA1 SNV gets an intronic cDNA notation."""
        # NC_000017.11:g.43045803C>A  (pipeline DB VRS: ga4gh:VA.lRBWfd_...)
        row = self._run('NC_000017.11:g.43045803C>A', self.BRCA1_MANE, hgvs_tools)
        assert row[HGVS_CDNA_COL] != '-'
        assert row[REFERENCE_SEQUENCE_COL] == self.BRCA1_MANE
        assert row[HGVS_CDNA_COL] == 'c.5468-1G>T'

    def test_brca1_intronic_deletion(self, hgvs_tools):
        """gnomAD-only intronic BRCA1 deletion gets an intronic cDNA notation."""
        # NC_000017.11:g.43045810del (approximated from pipeline row CAGAG>C)
        row = self._run('NC_000017.11:g.43045810del', self.BRCA1_MANE, hgvs_tools)
        assert row[HGVS_CDNA_COL] != '-'
        assert re.match(r'c\.\d+[+-]\d+', row[HGVS_CDNA_COL]), (
            f"Expected intronic cDNA notation, got: {row[HGVS_CDNA_COL]}")

    def test_brca1_near_5prime_utr(self, hgvs_tools):
        """gnomAD-only BRCA1 variant near 5' UTR gets a cDNA notation."""
        # NC_000017.11:g.43123988C>T  (pipeline DB VRS: ga4gh:VA.1InLp2...)
        row = self._run('NC_000017.11:g.43123988C>T', self.BRCA1_MANE, hgvs_tools)
        assert row[HGVS_CDNA_COL] != '-'
        assert row[REFERENCE_SEQUENCE_COL] == self.BRCA1_MANE

    def test_brca1_exonic_snv(self, hgvs_tools):
        """Exonic BRCA1 SNV (control): same result as before the fix."""
        row = self._run('NC_000017.11:g.43094692G>A', self.BRCA1_MANE, hgvs_tools)
        assert row[HGVS_CDNA_COL] == 'c.839C>T'
        assert row[REFERENCE_SEQUENCE_COL] == self.BRCA1_MANE
