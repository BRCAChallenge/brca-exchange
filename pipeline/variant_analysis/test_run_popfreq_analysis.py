"""
Unit tests for run_popfreq_analysis.py.

Pure-logic tests only -- no database connection is needed since all the
evidence-code computation (analyze_one_dataset / _compute_evidence_code) and
its supporting helpers (coverage/LCR parsing) are plain functions. main()
itself (the click entrypoint that talks to PostgreSQL) is intentionally not
exercised here.

Run with:
    pytest variant_analysis/test_run_popfreq_analysis.py
"""

import math

import pandas as pd
import pytest

import run_popfreq_analysis as popfreq


# ---------------------------------------------------------------------------
# field_defined
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("value,expected", [
    ("-", False),
    ("", True),
    ("0.001", True),
    (None, True),  # only the literal string "-" counts as undefined
])
def test_field_defined(value, expected):
    assert popfreq.field_defined(value) == expected


# ---------------------------------------------------------------------------
# read_lcr / overlaps_lcr
# ---------------------------------------------------------------------------

@pytest.fixture
def lcr_bed(tmp_path):
    bed_path = tmp_path / "lcr.bed"
    bed_path.write_text(
        "track name=lcr\n"
        "chr13\t100\t200\n"
        "chr13\t500\t600\n"
        "chr17\t50\t80\n"
    )
    return str(bed_path)


def test_read_lcr_parses_and_strips_chr_prefix(lcr_bed):
    lcr = popfreq.read_lcr(lcr_bed)
    assert set(lcr.keys()) == {"13", "17"}
    assert lcr["13"] == [(100, 200), (500, 600)]
    assert lcr["17"] == [(50, 80)]


def test_read_lcr_skips_header_and_short_lines(tmp_path):
    bed_path = tmp_path / "lcr.bed"
    bed_path.write_text("browser position chr13:1-1000\n#comment\nchr13\t100\n")
    lcr = popfreq.read_lcr(str(bed_path))
    assert lcr == {}


@pytest.mark.parametrize("chrom,start,end,expected", [
    ("13", 150, 150, True),     # fully inside chr13:100-200 (0-based half-open)
    ("13", 90, 100, False),     # 1-based inclusive [90,100] -> 0-based [89,100); ends exactly at
                                # the LCR region's start (100), which is excluded (half-open) -> no overlap
    ("13", 300, 400, False),    # between the two chr13 regions
    ("13", 600, 700, True),     # 1-based start 600 -> 0-based 599, which is inside [500,600) (half-open)
    ("17", 60, 60, True),
    ("18", 60, 60, False),      # chromosome not present in LCR at all
])
def test_overlaps_lcr(lcr_bed, chrom, start, end, expected):
    lcr = popfreq.read_lcr(lcr_bed)
    assert popfreq.overlaps_lcr(chrom, start, end, lcr) == expected


def test_overlaps_lcr_with_chr_prefixed_query_still_matches_stripped_keys(lcr_bed):
    # overlaps_lcr does not itself strip 'chr' -- callers are expected to pass
    # bare chromosome numbers, matching how read_lcr stores its keys.
    lcr = popfreq.read_lcr(lcr_bed)
    assert popfreq.overlaps_lcr("chr13", 150, 150, lcr) is False
    assert popfreq.overlaps_lcr("13", 150, 150, lcr) is True


# ---------------------------------------------------------------------------
# estimate_coverage
# ---------------------------------------------------------------------------

def test_estimate_coverage_observable_uses_mean_of_positions():
    cov_data = {
        1: {
            100: {"mean": 20.0, "median": 18.0},
            101: {"mean": 30.0, "median": 32.0},
        }
    }
    observable, coverage = popfreq.estimate_coverage(100, 101, 1, cov_data)
    assert observable is True
    assert coverage == pytest.approx(25.0)


def test_estimate_coverage_uses_min_of_mean_and_median_when_requested():
    cov_data = {1: {100: {"mean": 40.0, "median": 10.0}}}
    observable, coverage = popfreq.estimate_coverage(100, 100, 1, cov_data, use_median=True)
    assert observable is True
    assert coverage == pytest.approx(10.0)


def test_estimate_coverage_not_observable_when_chrom_absent():
    cov_data = {1: {100: {"mean": 20.0}}}
    observable, coverage = popfreq.estimate_coverage(100, 100, 2, cov_data)
    assert observable is False
    assert coverage == 0


def test_estimate_coverage_not_observable_when_positions_absent():
    cov_data = {1: {999: {"mean": 20.0}}}
    observable, coverage = popfreq.estimate_coverage(100, 100, 1, cov_data)
    assert observable is False
    assert coverage == 0


# ---------------------------------------------------------------------------
# read_coverage
# ---------------------------------------------------------------------------

def test_read_coverage_roundtrips_parquet(tmp_path):
    df = pd.DataFrame({
        "chrom": [13, 13, 17],
        "pos": [100, 101, 50],
        "mean": [20.0, 25.0, 30.0],
    })
    path = tmp_path / "coverage.parquet"
    df.to_parquet(path)

    coverage = popfreq.read_coverage(str(path))

    assert coverage[13][100] == {"mean": 20.0}
    assert coverage[13][101] == {"mean": 25.0}
    assert coverage[17][50] == {"mean": 30.0}


def test_read_coverage_falls_back_to_weighted_mean_coverage_column(tmp_path):
    df = pd.DataFrame({"chrom": [13], "pos": [100], "weighted_mean_coverage": [15.0]})
    path = tmp_path / "coverage.parquet"
    df.to_parquet(path)

    coverage = popfreq.read_coverage(str(path))

    assert coverage[13][100] == {"mean": 15.0}


def test_read_coverage_includes_median_when_present(tmp_path):
    df = pd.DataFrame({"chrom": [13], "pos": [100], "mean": [15.0], "median": [14.0]})
    path = tmp_path / "coverage.parquet"
    df.to_parquet(path)

    coverage = popfreq.read_coverage(str(path))

    assert coverage[13][100] == {"mean": 15.0, "median": 14.0}


# ---------------------------------------------------------------------------
# is_sufficient
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("read_depth,rare_variant,expected", [
    (25, True, True),
    (24, True, False),
    (20, False, True),
    (19, False, False),
])
def test_is_sufficient(read_depth, rare_variant, expected):
    assert popfreq.is_sufficient(read_depth, is_flagged=False, rare_variant=rare_variant) == expected


# ---------------------------------------------------------------------------
# analyze_one_dataset
# ---------------------------------------------------------------------------

def _config(**overrides):
    return popfreq.PopfreqConfig(**overrides)


def test_analyze_one_dataset_filter_flag_fails_regardless_of_everything_else():
    code, msg = popfreq.analyze_one_dataset(
        "0.5", "1000", True, 100, True, 1,
        "13", 100, 100, None, "0.5", "afr", "-", "-",
        config=_config(), gnomad_label="gnomAD v4.1", debug=False,
    )
    assert code == popfreq.FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG


@pytest.mark.parametrize("faf,expected_code", [
    ("0.01", popfreq.BA1),               # > ba1 threshold (0.001)
    ("0.0005", popfreq.BS1),             # between bs1 (0.0001) and ba1
    ("0.00005", popfreq.BS1_SUPPORTING), # between bs1_supporting (0.00001) and bs1
    # Note: NO_CODE is unreachable from this branch under default config,
    # since bs1_supporting_faf_threshold == rare_variant_faf_threshold
    # (0.00001) -- any FAF that reaches the not-rare branch already exceeds
    # bs1_supporting_faf_threshold.
])
def test_analyze_one_dataset_common_variant_faf_thresholds(faf, expected_code):
    # Common-variant (not-rare) path requires read_depth >= 20 and
    # faf > rare_variant_faf_threshold (0.00001) to stay on this branch.
    code, msg = popfreq.analyze_one_dataset(
        faf, "-", True, 20, False, 1,
        "13", 100, 100, None, faf, "afr", "-", "-",
        config=_config(), gnomad_label="gnomAD v4.1", debug=False,
    )
    assert code == expected_code
    assert "gnomAD v4.1" in msg


def test_analyze_one_dataset_faf_at_or_below_rare_threshold_is_rare_not_common():
    # A FAF at-or-below rare_variant_faf_threshold (0.00001) reclassifies the
    # variant as rare, which requires read_depth >= 25 (not 20) to proceed.
    faf = "0.000001"
    code, _ = popfreq.analyze_one_dataset(
        faf, "-", True, 20, False, 1,
        "13", 100, 100, None, faf, "afr", "-", "-",
        config=_config(), gnomad_label="gnomAD v4.1", debug=False,
    )
    assert code == popfreq.FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG

    code, _ = popfreq.analyze_one_dataset(
        faf, "-", True, 25, False, 1,
        "13", 100, 100, None, faf, "afr", "-", "-",
        config=_config(), gnomad_label="gnomAD v4.1", debug=False,
    )
    assert code == popfreq.PM2_SUPPORTING


def test_analyze_one_dataset_common_variant_insufficient_depth():
    code, msg = popfreq.analyze_one_dataset(
        "0.01", "-", True, 19, False, 1,
        "13", 100, 100, None, "0.01", "afr", "-", "-",
        config=_config(), gnomad_label="gnomAD v4.1", debug=False,
    )
    assert code == popfreq.FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG


def test_analyze_one_dataset_rare_variant_insufficient_depth():
    code, msg = popfreq.analyze_one_dataset(
        "-", "-", True, 24, False, 1,
        "13", 100, 100, None, "-", "-", "-", "-",
        config=_config(), gnomad_label="gnomAD v4.1", debug=False,
    )
    assert code == popfreq.FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG


def test_analyze_one_dataset_rare_with_faf_and_high_allele_count_is_no_code():
    # FAF *is* defined (just below the rare threshold) and allele count is
    # convincingly high -- this is not "no FAF computed", so it should always
    # be NO_CODE regardless of missing_faf_suggests_absence.
    for missing_faf_suggests_absence in (False, True):
        code, msg = popfreq.analyze_one_dataset(
            "0.000001", "1000", True, 25, False, 1,
            "13", 100, 100, None, "0.000001", "afr", "-", "-",
            config=_config(missing_faf_suggests_absence=missing_faf_suggests_absence),
            gnomad_label="gnomAD v4.1", debug=False,
        )
        assert code == popfreq.NO_CODE


def test_analyze_one_dataset_missing_allele_count_is_pm2_candidate():
    code, msg = popfreq.analyze_one_dataset(
        "-", "-", True, 25, False, 1,
        "13", 100, 100, None, "-", "-", "-", "-",
        config=_config(), gnomad_label="gnomAD v4.1", debug=False,
    )
    assert code == popfreq.PM2_SUPPORTING
    assert "gnomAD v4.1" in msg


def test_analyze_one_dataset_low_allele_count_is_pm2_candidate():
    code, msg = popfreq.analyze_one_dataset(
        "-", "1", True, 25, False, 1,  # allele_count_threshold=1, allele_count=1 -> not convincing
        "13", 100, 100, None, "-", "-", "-", "-",
        config=_config(), gnomad_label="gnomAD v4.1", debug=False,
    )
    assert code == popfreq.PM2_SUPPORTING


def test_no_code_message_quotes_the_faf_when_one_was_computed():
    code, msg = popfreq.analyze_one_dataset(
        "0.000001", "1000", True, 25, False, 1,
        "13", 100, 100, None, "0.000001", "afr", "-", "-",
        config=_config(), gnomad_label="gnomAD v4.1", debug=False,
    )
    assert code == popfreq.NO_CODE
    assert "is 0.000001 in the afr genetic ancestry group" in msg
    assert "not calculated" not in msg
    assert "95% CI" in msg and "%%" not in msg


def test_no_code_message_says_not_calculated_when_faf_is_missing():
    code, msg = popfreq.analyze_one_dataset(
        "-", "1000", True, 25, False, 1,
        "13", 100, 100, None, "-", "-", "-", "-",
        config=_config(), gnomad_label="gnomAD v4.1", debug=False,
    )
    assert code == popfreq.NO_CODE
    assert "was not calculated" in msg
    assert "is - in the -" not in msg


# --- missing_faf_suggests_absence: the popfreq_1.2 vs popfreq_1.3 behavior ---

def test_missing_faf_high_allele_count_blocked_by_default_legacy_behavior():
    """popfreq_1.2 (default config): missing FAF + high allele count -> NO_CODE."""
    code, msg = popfreq.analyze_one_dataset(
        "-", "1000", True, 25, False, 1,
        "13", 100, 100, None, "-", "-", "-", "-",
        config=_config(), gnomad_label="gnomAD v4.1", debug=False,
    )
    assert code == popfreq.NO_CODE


def test_missing_faf_high_allele_count_is_pm2_candidate_under_1_3():
    """popfreq_1.3 (missing_faf_suggests_absence=True): missing FAF alone makes
    the variant a PM2_Supporting candidate no matter how high the allele count is."""
    code, msg = popfreq.analyze_one_dataset(
        "-", "1000000", True, 25, False, 1,
        "13", 100, 100, None, "-", "-", "-", "-",
        config=_config(missing_faf_suggests_absence=True),
        gnomad_label="gnomAD v4.1", debug=False,
    )
    assert code == popfreq.PM2_SUPPORTING


def test_missing_faf_bypass_still_respects_indel_size_gate():
    """The 1.3 bypass makes the variant a *candidate*; a large indel still
    fails the snv_or_small_indel eligibility check and gets NO_CODE_INDEL."""
    code, msg = popfreq.analyze_one_dataset(
        "-", "1000000", False, 25, False, 1,  # snv_or_small_indel=False
        "13", 100, 100, None, "-", "-", "-", "-",
        config=_config(missing_faf_suggests_absence=True),
        gnomad_label="gnomAD v4.1", debug=False,
    )
    assert code == popfreq.NO_CODE_INDEL


def test_missing_faf_bypass_still_respects_lcr_gate():
    """The 1.3 bypass makes the variant a *candidate*; overlapping an LCR
    still fails and gets FAIL_LCR."""
    lcr = {"13": [(99, 101)]}
    code, msg = popfreq.analyze_one_dataset(
        "-", "1000000", True, 25, False, 1,
        "13", 100, 100, lcr, "-", "-", "-", "-",
        config=_config(missing_faf_suggests_absence=True, use_lcr=True),
        gnomad_label="gnomAD v4.1", debug=False,
    )
    assert code == popfreq.FAIL_LCR


def test_missing_allele_count_suggests_absence_can_be_disabled():
    """Setting missing_allele_count_suggests_absence=False makes a missing
    allele count no longer imply absence (mirrors missing_faf_suggests_absence's
    default-off semantics, exercised for symmetry/coverage of the setting)."""
    code, msg = popfreq.analyze_one_dataset(
        "0.000001", "-", True, 25, False, 1,
        "13", 100, 100, None, "0.000001", "afr", "-", "-",
        config=_config(missing_allele_count_suggests_absence=False),
        gnomad_label="gnomAD v4.1", debug=False,
    )
    # missing allele count no longer forces "not convincingly present" on its
    # own, but FAF *is* defined here, so the field_defined(faf) branch is hit.
    assert code == popfreq.NO_CODE


# ---------------------------------------------------------------------------
# _compute_evidence_code: coverage-dataset fallback + gnomAD-version naming
# ---------------------------------------------------------------------------

_COORDS = ("hgvs.test", "chr13", 100, "A", "G")

# Frequency rows are keyed the way report_gnomad labels releases -- note 'v3',
# not 'v3.1', which is the mapping _COVERAGE_TO_REPORT exists to bridge.
_JOINT = ("v4.1", "joint")
_GENOME_REPORT = ("v3", "genome")


def _freq(flags=None, allele_count="5", faf95="0.01", faf95_pop="afr"):
    return (flags, allele_count, faf95, faf95_pop)


def test_compute_evidence_code_prefers_first_sufficient_dataset():
    coverage_v4_joint = {13: {100: {"mean": 30.0}}}  # sufficient
    coverage_v3_genome = {13: {100: {"mean": 30.0}}}  # also sufficient, but lower priority
    coverage_meta = [("v4.1", "joint"), ("v3.1", "genome")]

    code, msg, idx = popfreq._compute_evidence_code(
        *_COORDS, {_JOINT: _freq(), _GENOME_REPORT: _freq()},
        [coverage_v4_joint, coverage_v3_genome], coverage_meta,
        None, _config(),
    )
    assert idx == 0
    assert "gnomAD v4.1" in msg
    assert "gnomAD v3.1" not in msg


def test_compute_evidence_code_falls_back_to_next_dataset_when_first_insufficient():
    coverage_v4_joint = {13: {100: {"mean": 2.0}}}   # insufficient
    coverage_v3_genome = {13: {100: {"mean": 30.0}}}  # sufficient
    coverage_meta = [("v4.1", "joint"), ("v3.1", "genome")]

    code, msg, idx = popfreq._compute_evidence_code(
        *_COORDS, {_JOINT: _freq(), _GENOME_REPORT: _freq()},
        [coverage_v4_joint, coverage_v3_genome], coverage_meta,
        None, _config(),
    )
    assert idx == 1
    assert "gnomAD v3.1" in msg
    assert "gnomAD v4.1" not in msg


def test_compute_evidence_code_no_sufficient_dataset_fails_with_no_dataset_idx():
    coverage_v4_joint = {13: {100: {"mean": 1.0}}}
    coverage_meta = [("v4.1", "joint")]

    code, msg, idx = popfreq._compute_evidence_code(
        *_COORDS, {_JOINT: _freq()}, [coverage_v4_joint], coverage_meta, None, _config(),
    )
    assert code == popfreq.FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG
    assert idx is None


def test_compute_evidence_code_row_with_no_frequency_data_at_all_is_skipped():
    """A release whose row carries neither an allele count nor a FAF is a
    loading gap, so it is skipped rather than scored as an absence. Absence is
    represented by having no row for the release at all (see
    test_release_with_no_report_row_but_good_coverage_is_a_genuine_absence).

    When the skipped release is the only candidate there is nothing left to
    score against, so the read-depth failure code is what comes back.
    """
    coverage_v4_joint = {13: {100: {"mean": 30.0}}}
    coverage_meta = [("v4.1", "joint")]

    code, msg, idx = popfreq._compute_evidence_code(
        *_COORDS,
        {_JOINT: _freq(flags=None, allele_count=None, faf95=None, faf95_pop=None)},
        [coverage_v4_joint], coverage_meta, None, _config(),
    )
    assert code == popfreq.FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG
    assert idx is None


def test_compute_evidence_code_end_to_end_missing_faf_suggests_absence():
    """Exercises the full popfreq_1.3 request end-to-end via _compute_evidence_code:
    a variant absent a FAF value but with a large recorded allele count is a
    PM2_Supporting candidate when missing_faf_suggests_absence=True."""
    coverage_v4_joint = {13: {100: {"mean": 30.0}}}
    coverage_meta = [("v4.1", "joint")]

    code_legacy, _, _ = popfreq._compute_evidence_code(
        *_COORDS, {_JOINT: _freq(allele_count="9999", faf95=None)},
        [coverage_v4_joint], coverage_meta, None, _config(),
    )
    assert code_legacy == popfreq.NO_CODE

    code_1_3, _, _ = popfreq._compute_evidence_code(
        *_COORDS, {_JOINT: _freq(allele_count="9999", faf95=None)},
        [coverage_v4_joint], coverage_meta, None,
        _config(missing_faf_suggests_absence=True),
    )
    assert code_1_3 == popfreq.PM2_SUPPORTING


# ---------------------------------------------------------------------------
# The selected release supplies coverage AND frequencies (the popfreq bug fix)
# ---------------------------------------------------------------------------

def test_frequencies_come_from_the_release_that_supplied_coverage():
    """The whole point of the fix: when v4.1 joint coverage is insufficient and
    the v3.1 genome release is used instead, the FAF must be read from the v3
    report row -- not from v4.1 joint, which used to be hardcoded."""
    coverage_v4_joint = {13: {100: {"mean": 2.0}}}    # insufficient
    coverage_v3_genome = {13: {100: {"mean": 30.0}}}  # sufficient
    coverage_meta = [("v4.1", "joint"), ("v3.1", "genome")]

    code, msg, idx = popfreq._compute_evidence_code(
        *_COORDS,
        {_JOINT: _freq(faf95="0.01"),           # would be BA1 if wrongly used
         _GENOME_REPORT: _freq(faf95="0.0005")},  # BS1 band -- the correct source
        [coverage_v4_joint, coverage_v3_genome], coverage_meta, None, _config(),
    )
    assert idx == 1
    assert code == popfreq.BS1
    assert "0.0005" in msg and "gnomAD v3.1" in msg


def test_v3_1_coverage_label_maps_to_v3_report_rows():
    """coverage_meta says 'v3.1'; report_gnomad says 'v3'. Frequencies keyed the
    report way must still be found."""
    coverage = {13: {100: {"mean": 30.0}}}
    code, msg, idx = popfreq._compute_evidence_code(
        *_COORDS, {_GENOME_REPORT: _freq(faf95="0.01")},
        [coverage], [("v3.1", "genome")], None, _config(),
    )
    assert idx == 0
    assert code == popfreq.BA1   # found the row; not treated as absent


def test_rarity_is_evaluated_per_release_not_once_up_front():
    """Rarity selects the read-depth threshold (25 rare vs 20 common), so it has
    to be decided from each release's own FAF. Here v4.1 joint is common at
    depth 20 (sufficient), while judging it by the v3 release's rare FAF would
    have demanded 25 and wrongly skipped it."""
    coverage_v4_joint = {13: {100: {"mean": 20.0}}}
    coverage_v3_genome = {13: {100: {"mean": 30.0}}}
    coverage_meta = [("v4.1", "joint"), ("v3.1", "genome")]

    code, msg, idx = popfreq._compute_evidence_code(
        *_COORDS,
        {_JOINT: _freq(faf95="0.01"),               # common -> threshold 20 -> OK
         _GENOME_REPORT: _freq(faf95="0.000001")},  # rare
        [coverage_v4_joint, coverage_v3_genome], coverage_meta, None, _config(),
    )
    assert idx == 0
    assert code == popfreq.BA1


def test_release_with_no_frequency_data_loaded_is_skipped_not_read_as_absent():
    """The v4.1 exome rows exist but carry no AC and no FAF. Scoring them as an
    absence would manufacture PM2_Supporting, so the release is skipped and the
    next one is used."""
    coverage_exome = {13: {100: {"mean": 30.0}}}
    coverage_v3_genome = {13: {100: {"mean": 30.0}}}
    coverage_meta = [("v4.1", "exome"), ("v3.1", "genome")]

    code, msg, idx = popfreq._compute_evidence_code(
        *_COORDS,
        {("v4.1", "exome"): _freq(flags="-", allele_count="-", faf95="-", faf95_pop="-"),
         _GENOME_REPORT: _freq(faf95="0.01")},
        [coverage_exome, coverage_v3_genome], coverage_meta, None, _config(),
    )
    assert idx == 1
    assert code == popfreq.BA1


def test_release_with_no_report_row_but_good_coverage_is_a_genuine_absence():
    """No row at all means the variant was not observed in that release. With
    sufficient coverage that is a real absence and stays PM2-eligible."""
    coverage = {13: {100: {"mean": 30.0}}}
    code, msg, idx = popfreq._compute_evidence_code(
        *_COORDS, {}, [coverage], [("v4.1", "joint")], None, _config(),
    )
    assert idx == 0
    assert code == popfreq.PM2_SUPPORTING


# ---------------------------------------------------------------------------
# PopfreqConfig defaults -- lock in the documented popfreq_1.2/1.3 semantics
# ---------------------------------------------------------------------------

def test_popfreq_config_defaults_match_legacy_popfreq_1_2_behavior():
    config = popfreq.PopfreqConfig()
    assert config.missing_faf_suggests_absence is False
    assert config.missing_allele_count_suggests_absence is True
    assert config.use_lcr is True
    assert config.ba1_faf_threshold == popfreq._BA1_FAF_THRESHOLD
    assert config.bs1_faf_threshold == popfreq._BS1_FAF_THRESHOLD
