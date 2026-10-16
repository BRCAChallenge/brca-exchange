#!/usr/bin/env python
# coding: utf-8

"""
Populate analysis_provisional_evidence_codes from the database.

Reads gnomAD data from report_gnomad, GRCh38 coordinates from
variant_genomic_coordinates, and computes provisional population frequency
evidence codes for each variant. Results are upserted into
analysis_provisional_evidence_codes.

Releases are assessed in priority order: gnomAD v4.1 joint, then gnomAD v4.1
exome, then gnomAD v3.1 genome. The first release whose coverage is sufficient
for the variant is used, and that same release supplies the frequencies (FAF,
allele count, filter flags) the evidence code is computed from -- coverage and
frequency always come from the same gnomAD release. If no release is
sufficient, the evidence code is set to
FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG.

An LCR BED file may optionally be supplied as a local file.
"""

import bisect
import dataclasses
import logging
import math
import os
from typing import Optional

import click
import psycopg2
import psycopg2.extras

BA1 = "BA1 (met)"
BS1 = "BS1 (met)"
BS1_SUPPORTING = "BS1_Supporting (met)"
NO_CODE = "No code met (below threshold)"
NO_CODE_INDEL = "No code met (indel)"
PM2_SUPPORTING = "PM2_Supporting (met)"
FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG = "No code met (read depth, flags)"
FAIL_LCR = "No code met (low-complexity region)"

READ_DEPTH_THRESHOLD_FREQUENT_VARIANT = 20
READ_DEPTH_THRESHOLD_RARE_VARIANT = 25
ALLELE_COUNT_RARE_VARIANT_THRESHOLD = 1

# Coverage datasets are named after the gnomAD release ('v3.1' for the v3.1
# genome release); report_gnomad stores that same release under version 'v3'.
# The coverage-side name is what gets reported and stored, so this map is used
# only to find the release's frequency row.
_COVERAGE_TO_REPORT = {
    ('v4.1', 'joint'):  ('v4.1', 'joint'),
    ('v4.1', 'exome'):  ('v4.1', 'exome'),
    ('v3.1', 'genome'): ('v3',   'genome'),
}

# Default threshold values (popfreq_1.3)
_BA1_FAF_THRESHOLD            = 0.001
_BS1_FAF_THRESHOLD            = 0.0001
_BS1_SUPPORTING_FAF_THRESHOLD = 0.00001
_RARE_VARIANT_FAF_THRESHOLD   = 0.00001
_SMALL_INDEL_SIZE_THRESHOLD   = 50


@dataclasses.dataclass
class PopfreqConfig:
    """Configurable thresholds for population frequency evidence code computation."""
    ba1_faf_threshold:            float = _BA1_FAF_THRESHOLD
    bs1_faf_threshold:            float = _BS1_FAF_THRESHOLD
    bs1_supporting_faf_threshold: float = _BS1_SUPPORTING_FAF_THRESHOLD
    rare_variant_faf_threshold:   float = _RARE_VARIANT_FAF_THRESHOLD
    small_indel_size_threshold:   int   = _SMALL_INDEL_SIZE_THRESHOLD
    allele_count_threshold:       int   = ALLELE_COUNT_RARE_VARIANT_THRESHOLD
    use_lcr:                      bool  = True
    # Whether a missing allele count (variant not observed at all) suggests the
    # variant is absent from gnomAD, i.e. not convincingly present. Unchanged
    # from historical behavior; exposed here as an explicit, named setting.
    missing_allele_count_suggests_absence: bool = True
    # Whether a variant with no computed FAF at all is treated as suggesting
    # absence, so it remains a PM2_Supporting candidate regardless of allele
    # count. False (the legacy popfreq_1.2 behavior) by default -- an allele
    # count above allele_count_threshold blocks PM2_Supporting even when no
    # FAF was computed; popfreq_1.3 opts in explicitly via --missing-faf-suggests-absence.
    missing_faf_suggests_absence: bool = False

    def ba1_msg(self, faf, pop, ac, an, gnomad_label):
        return (
            f"The Total GrpMax filtering allele frequency (the lower threshold of the 95% CI) "
            f"in {gnomad_label} is {faf} in the {pop} genetic ancestry group (based on {ac}/{an} alleles) "
            f"which is above the ENIGMA BRCA1/2 VCEP threshold (>{self.ba1_faf_threshold}) for BA1 (BA1 met)."
        )

    def bs1_msg(self, faf, pop, ac, an, gnomad_label):
        return (
            f"The Total GrpMax filtering allele frequency (the lower threshold of the 95% CI) "
            f"in {gnomad_label} is {faf} in the {pop} genetic ancestry group (based on {ac}/{an} alleles) "
            f"which is above the ENIGMA BRCA1/2 VCEP threshold (>{self.bs1_faf_threshold}) for BS1, "
            f"and below the BA1 threshold (>{self.ba1_faf_threshold}) (BS1 met)."
        )

    def bs1_supporting_msg(self, faf, pop, ac, an, gnomad_label):
        return (
            f"The Total GrpMax filtering allele frequency (the lower threshold of the 95% CI) "
            f"in {gnomad_label} is {faf} in the {pop} genetic ancestry group (based on {ac}/{an} alleles) "
            f"which is above the ENIGMA BRCA1/2 VCEP threshold (>{self.rare_variant_faf_threshold}) "
            f"for BS1_Supporting, and below the BS1 threshold (>{self.bs1_faf_threshold}) (BS1_Supporting met)."
        )

    def no_code_met_msg(self, faf, pop, ac, an, gnomad_label):
        return (
            f"The Total GrpMax filtering allele frequency (the lower threshold of the 95% CI) "
            f"in {gnomad_label} is {faf} in the {pop} genetic ancestry group (based on {ac}/{an} alleles) "
            f"which is below the ENIGMA BRCA1/2 VCEP threshold (>{self.bs1_supporting_faf_threshold}) "
            f"for BS1_Supporting and does not meet any population code "
            f"(BA1, BS1, BS1_Supporting, PM2_Supporting are not met)."
        )


def no_code_no_faf_msg(gnomad_label):
    return (
        f"This variant is recorded in {gnomad_label}, however the Total GrpMax filtering allele frequency "
        "(the lower threshold of the 95% CI) was not calculated, therefore this variant does not meet "
        "any population code (BA1, BS1, BS1_Supporting, PM2_Supporting are not met)."
    )


def no_code_indel_msg(gnomad_label):
    return (
        "This [duplication/insertion/deletion/delins/large genomic rearrangement] variant "
        f"was not observed in {gnomad_label}, but PM2_Supporting was not applied since recall "
        "is considered suboptimal for this type of variant (PM2_Supporting not met). "
    )


def pm2_supporting_msg(gnomad_label):
    return f"This variant is absent from {gnomad_label} (PM2_Supporting met)."


FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG_MSG = (
    f"This variant is present in gnomAD but is not meeting the "
    f"specified read depths threshold ≥{READ_DEPTH_THRESHOLD_FREQUENT_VARIANT} "
    f"(PM2_Supporting, BS1, and BA1 are not met)."
)
FAIL_LCR_MSG = "PM2_Supporting was not applied since this variant overlaps a low-complexity region (LCR)."

DB_BATCH = 500
log = logging.getLogger(__name__)


def read_lcr(lcr_file):
    """
    Read a BED file of low-complexity regions.
    Returns: dict with structure {chrom: sorted list of (start, end) tuples}
    Coordinates are 0-based half-open, as per BED format.
    'chr' prefixes are stripped from chromosome names.
    """
    lcr = {}
    with open(lcr_file) as f:
        for line in f:
            if line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            chrom = fields[0].lstrip('chr')
            start, end = int(fields[1]), int(fields[2])
            lcr.setdefault(chrom, []).append((start, end))
    for chrom in lcr:
        lcr[chrom].sort()
    return lcr


def overlaps_lcr(chrom, start, end, lcr):
    """
    Return True if the variant region overlaps any low-complexity region.
    Variant coords (start, end) are 1-based inclusive (VCF-style).
    LCR regions are 0-based half-open (BED-style).
    """
    chrom_key = str(chrom)
    if chrom_key not in lcr:
        return False
    regions = lcr[chrom_key]
    var_start = start - 1
    var_end = end
    starts = [r[0] for r in regions]
    idx = bisect.bisect_left(starts, var_end)
    for i in range(idx - 1, -1, -1):
        r_start, r_end = regions[i]
        if r_end <= var_start:
            break
        return True
    return False


def read_coverage(coverage_parquet):
    """
    Read coverage data from a Parquet file and organize by chromosome and position.
    Returns: dict with structure {chrom: {pos: {"mean": float[, "median": float]}}}
    Expected columns: chrom (int), pos (int), and one of: mean, weighted_mean_coverage.
    Optional: median.
    """
    import pandas as pd
    df = pd.read_parquet(coverage_parquet)
    mean_col = 'mean' if 'mean' in df.columns else 'weighted_mean_coverage'
    has_median = 'median' in df.columns
    coverage = {}
    for chrom_val, grp in df.groupby(df['chrom'].astype(int)):
        pos_arr = grp['pos'].astype(int).values
        mean_arr = grp[mean_col].astype(float).values
        if has_median:
            med_arr = grp['median'].astype(float).values
            coverage[int(chrom_val)] = {
                int(p): {'mean': float(m), 'median': float(md)}
                for p, m, md in zip(pos_arr, mean_arr, med_arr)
            }
        else:
            coverage[int(chrom_val)] = {
                int(p): {'mean': float(m)}
                for p, m in zip(pos_arr, mean_arr)
            }
    return coverage


def estimate_coverage(start, end, chrom, cov_data, debug=False, use_median=False):
    """
    Estimate coverage for a genomic region from coverage data stored as a dictionary.
    cov_data: dict with structure {chrom: {pos: {"mean": float, "median": float}}}
    """
    positions = list(range(start, end+1))
    chrom_key = int(chrom)

    mean_values = []
    median_values = []

    if chrom_key in cov_data:
        for pos in positions:
            if pos in cov_data[chrom_key]:
                mean_values.append(cov_data[chrom_key][pos]["mean"])
                if use_median:
                    median_values.append(cov_data[chrom_key][pos]["median"])

    if len(mean_values) > 0:
        observable = True
        meanval = sum(mean_values) / len(mean_values)
        if not use_median:
            coverage = meanval
            medianval = None
        else:
            median_values_sorted = sorted(median_values)
            n = len(median_values_sorted)
            if n % 2 == 0:
                medianval = (median_values_sorted[n//2-1] + median_values_sorted[n//2]) / 2
            else:
                medianval = median_values_sorted[n//2]
            coverage = min(meanval, medianval)
    else:
        observable = False
        coverage = 0
        meanval = None
        medianval = None

    if debug:
        print("coverage assessment: observable:", observable, "meanval:",
              meanval, "coverage", coverage)
    return(observable, coverage)


def field_defined(field):
    return field != "-"


def is_sufficient(read_depth, is_flagged, rare_variant):
    """
    Return True if a coverage dataset is sufficient for a variant.

    Sufficiency is determined by read depth alone — the flag check is handled
    downstream in analyze_one_dataset, which assigns the appropriate code while
    still recording which coverage dataset was used.
    """
    threshold = (READ_DEPTH_THRESHOLD_RARE_VARIANT if rare_variant
                 else READ_DEPTH_THRESHOLD_FREQUENT_VARIANT)
    return read_depth >= threshold


def analyze_one_dataset(faf95_popmax_str, allele_count, snv_or_small_indel,
                        read_depth, vcf_filter_flag,
                        allele_count_threshold,
                        chrom, genome_start, genome_end, lcr,
                        faf95_popmax, faf95_popmax_population,
                        allele_count_pop, allele_number_pop,
                        config: PopfreqConfig, gnomad_label, debug=True):
    if vcf_filter_flag:
        return(FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG, FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG_MSG)
    rare_variant = False
    if field_defined(faf95_popmax_str):
        faf = float(faf95_popmax_str)
        if math.isnan(faf):
            rare_variant = True
        elif faf <= config.rare_variant_faf_threshold:
            rare_variant = True
    else:
        faf = None
        rare_variant = True
    if debug:
        print("Rare variant", rare_variant, "read depth", read_depth, "flagged", vcf_filter_flag)
    if rare_variant and read_depth < READ_DEPTH_THRESHOLD_RARE_VARIANT:
        return(FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG,
               FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG_MSG)
    if (not rare_variant) and read_depth < READ_DEPTH_THRESHOLD_FREQUENT_VARIANT:
        return(FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG,
               FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG_MSG)
    if not rare_variant:
        if debug:
            print("Not rare variant.  FAF:", faf)
        if faf > config.ba1_faf_threshold:
            return(BA1, config.ba1_msg(faf95_popmax, faf95_popmax_population, allele_count_pop, allele_number_pop, gnomad_label))
        elif faf > config.bs1_faf_threshold:
            return(BS1, config.bs1_msg(faf95_popmax, faf95_popmax_population, allele_count_pop, allele_number_pop, gnomad_label))
        elif faf > config.bs1_supporting_faf_threshold:
            return(BS1_SUPPORTING, config.bs1_supporting_msg(faf95_popmax, faf95_popmax_population, allele_count_pop, allele_number_pop, gnomad_label))
        else:
            return(NO_CODE, config.no_code_met_msg(faf95_popmax, faf95_popmax_population, allele_count_pop, allele_number_pop, gnomad_label))
    if debug:
        print("Rare variant.  Allele count", allele_count, "SNV", snv_or_small_indel)
    if not field_defined(allele_count):
        variant_convincingly_present = not config.missing_allele_count_suggests_absence
    elif int(allele_count) <= allele_count_threshold:
        variant_convincingly_present = False
    else:
        variant_convincingly_present = True
    if not field_defined(faf95_popmax_str) and config.missing_faf_suggests_absence:
        variant_convincingly_present = False
    if variant_convincingly_present:
        if field_defined(faf95_popmax_str):
            return(NO_CODE, config.no_code_met_msg(faf95_popmax, faf95_popmax_population, allele_count_pop, allele_number_pop, gnomad_label))
        else:
            return(NO_CODE, no_code_no_faf_msg(gnomad_label))
    if config.use_lcr and lcr is not None:
        if overlaps_lcr(chrom, genome_start, genome_end, lcr):
            return(FAIL_LCR, FAIL_LCR_MSG)
    if snv_or_small_indel:
        return(PM2_SUPPORTING, pm2_supporting_msg(gnomad_label))
    else:
        return(NO_CODE_INDEL, no_code_indel_msg(gnomad_label))


def _as_field(value):
    """NULL gnomAD columns (variant not reported in that release) read as '-'."""
    return value if value is not None else '-'


def _compute_evidence_code(hgvs_cdna, chr_, pos, ref, alt,
                            freq_by_dataset,
                            coverage_datasets, coverage_meta, lcr, config: PopfreqConfig, debug=False):
    """Compute the popfreq evidence code for one variant from DB row fields.

    coverage_datasets is a list of coverage dicts tried in priority order.
    coverage_meta is the parallel list of (gnomad_version, gnomad_data_type)
    tuples naming each release.

    freq_by_dataset maps (version, data_type) -- keyed as report_gnomad labels
    them, see _COVERAGE_TO_REPORT -- to that release's
    (Flags, Allele_count, faf95_popmax, faf95_popmax_population).

    The release that supplies the read depth also supplies the frequencies:
    each candidate's own FAF decides whether the variant is rare there (which
    selects the read-depth threshold), and its own coverage decides whether it
    is usable at all. The first usable release wins. If none is usable,
    FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG is returned.
    """
    start = int(pos)
    end = start + len(ref) - 1
    chrom = int(chr_.lstrip('chr'))

    is_snv = (len(ref) == 1 and len(alt) == 1)
    indel_size = max(len(ref), len(alt)) - 1
    snv_or_small_indel = is_snv or (indel_size <= config.small_indel_size_threshold)

    if debug:
        print(f'  _compute_evidence_code inputs:')
        print(f'    hgvs_cdna={hgvs_cdna!r}  chr={chr_}  pos={pos}  ref={ref!r}')
        print(f'    frequencies available for: {sorted(freq_by_dataset)}')
        print(f'  derived:')
        print(f'    start={start}  end={end}  chrom={chrom}  is_snv={is_snv}  indel_size={indel_size}  snv_or_small_indel={snv_or_small_indel}')

    for idx, cov_data in enumerate(coverage_datasets):
        row = freq_by_dataset.get(_COVERAGE_TO_REPORT[coverage_meta[idx]])

        if row is None:
            # Not reported in this release. With sufficient coverage that is a
            # genuine absence, so the release stays in the running.
            flags, allele_count, faf95, faf95_pop = None, None, None, None
        else:
            flags, allele_count, faf95, faf95_pop = row
            if not field_defined(_as_field(allele_count)) and not field_defined(_as_field(faf95)):
                # The row exists but carries no frequency data at all. That is a
                # gap in what was loaded, not an observed absence -- scoring it
                # as absent would manufacture PM2_Supporting -- so skip the
                # release entirely and fall through to the next one.
                log.warning('%s: no frequency data loaded for gnomAD %s %s; skipping that release',
                            hgvs_cdna, *coverage_meta[idx])
                continue

        faf95 = _as_field(faf95)
        allele_count = _as_field(allele_count)
        faf95_pop = _as_field(faf95_pop)
        # report_gnomad.Flags uses '-' for PASS; any other non-null value is a filter flag.
        is_flagged = flags is not None and flags != '-'

        # Rarity is decided by this release's own FAF (mirrors analyze_one_dataset),
        # and sets which read-depth threshold applies.
        if field_defined(faf95):
            faf_val = float(faf95)
            rare_variant = math.isnan(faf_val) or faf_val <= config.rare_variant_faf_threshold
        else:
            rare_variant = True

        _, depth = estimate_coverage(start, end, chrom, cov_data, debug=debug)

        if debug:
            print(f'    candidate {coverage_meta[idx]}: faf95={faf95!r} allele_count={allele_count!r} '
                  f'faf95_pop={faf95_pop!r} is_flagged={is_flagged} rare_variant={rare_variant} depth={depth}')

        if not is_sufficient(depth, is_flagged, rare_variant):
            continue

        gnomad_version, _gnomad_data_type = coverage_meta[idx]
        gnomad_label = f"gnomAD {gnomad_version}"

        # Per-population allele counts are not stored in the DB; '-' is used as placeholder
        # in the message text only (does not affect code assignment logic).
        code, msg = analyze_one_dataset(
            faf95, allele_count, snv_or_small_indel, depth, is_flagged,
            config.allele_count_threshold,
            chrom, start, end, lcr,
            faf95, faf95_pop,
            '-', '-',
            config=config,
            gnomad_label=gnomad_label,
            debug=debug,
        )

        if debug:
            print(f'  output: code={code!r} from {coverage_meta[idx]}')

        return code, msg, idx

    if debug:
        print(f'  output: code={FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG!r} (no sufficient dataset)')
    return FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG, FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG_MSG, None


_QUERY_ALL = """
SELECT v."VRS_Digest", v."HGVS_cDNA",
       gc.chr, gc.pos, gc.ref, gc.alt
FROM variant v
JOIN variant_genomic_coordinates gc
     ON gc."VRS_Digest" = v."VRS_Digest" AND gc.assembly = 'GRCh38'
"""

# Unscored is per method: a variant already carrying, say, popfreq_1.3 is still
# unscored for popfreq_1.2.
_QUERY_UNSCORED = _QUERY_ALL + """
LEFT JOIN analysis_provisional_evidence_codes apec
     ON apec."VRS_Digest" = v."VRS_Digest"
    AND apec.method_name IS NOT DISTINCT FROM %s
WHERE apec."VRS_Digest" IS NULL
"""

_QUERY_ONE = _QUERY_ALL + """
WHERE v."VRS_Digest" = %s
"""

# Every release's frequencies, keyed by (version, data_type) as report_gnomad
# labels them. The release selected for coverage supplies the frequencies too.
_QUERY_FREQ = """
SELECT "VRS_Digest", version, data_type,
       "Flags", "Allele_count", faf95_popmax, faf95_popmax_population
FROM report_gnomad
"""

_QUERY_FREQ_ONE = _QUERY_FREQ + """
WHERE "VRS_Digest" = %s
"""


@click.command()
@click.option('--db-url', default='postgresql://postgres:postgres@localhost/storage.pg',
              envvar='PIPELINE_DB_URL', show_default=True)
@click.option('--schema', default='pipeline', show_default=True)
@click.option('--coverage-v4-joint', required=True, type=click.Path(exists=True),
              help='Parquet file of gnomAD v4.1 joint (exome+genome) per-position coverage')
@click.option('--coverage-v4-exome', default=None, type=click.Path(exists=True),
              help='Parquet file of gnomAD v4.1 exome per-position coverage (fallback)')
@click.option('--coverage-v3-genome', default=None, type=click.Path(exists=True),
              help='Parquet file of gnomAD v3.1 genome per-position coverage (fallback)')
@click.option('--lcr', default=None, type=click.Path(exists=True),
              help='BED file of low-complexity regions; variants overlapping an LCR '
                   'will not be assigned PM2_Supporting')
@click.option('--overwrite', is_flag=True, default=False,
              help='Re-score variants already in analysis_provisional_evidence_codes')
@click.option('--debug', is_flag=True, default=False,
              help='Print per-variant debugging info; implies --dry-run')
@click.option('--dry-run', is_flag=True, default=False,
              help='Compute and print results without writing to the database')
@click.option('--vrs-digest', default=None, metavar='DIGEST',
              help='Analyze only the single variant with this VRS digest')
@click.option('--method-name', default=None,
              help='Method name to record in the method_name column (e.g. "popfreq_1.3")')
@click.option('--bs1-supporting-faf-threshold', default=_BS1_SUPPORTING_FAF_THRESHOLD, type=float,
              show_default=True, help='FAF threshold above which BS1_Supporting is assigned')
@click.option('--rare-variant-faf-threshold', default=_RARE_VARIANT_FAF_THRESHOLD, type=float,
              show_default=True, help='FAF threshold at-or-below which a variant is considered rare')
@click.option('--small-indel-size-threshold', default=_SMALL_INDEL_SIZE_THRESHOLD, type=int,
              show_default=True, help='Max ref length (bp) to be treated as a small indel eligible for PM2_Supporting')
@click.option('--allele-count-rare-variant-threshold', default=ALLELE_COUNT_RARE_VARIANT_THRESHOLD, type=int,
              show_default=True, help='Allele count above which PM2_Supporting is not assigned for rare variants')
@click.option('--no-lcr', is_flag=True, default=False,
              help='Disable LCR filtering (PM2_Supporting may be assigned regardless of LCR overlap)')
@click.option('--missing-faf-suggests-absence', is_flag=True, default=False,
              help='popfreq_1.3 behavior: treat a variant with no computed FAF at all as suggesting '
                   'absence, so it remains a PM2_Supporting candidate regardless of allele count. '
                   'By default (legacy popfreq_1.2 behavior), an allele count above the rare-variant '
                   'threshold blocks PM2_Supporting even when no FAF was computed.')
def main(db_url, schema, coverage_v4_joint, coverage_v4_exome, coverage_v3_genome,
         lcr, overwrite, debug, dry_run, vrs_digest, method_name,
         bs1_supporting_faf_threshold, rare_variant_faf_threshold,
         small_indel_size_threshold, allele_count_rare_variant_threshold, no_lcr,
         missing_faf_suggests_absence):
    logging.basicConfig(level=logging.WARNING, format='%(levelname)s %(message)s')

    if debug:
        dry_run = True

    config = PopfreqConfig(
        bs1_supporting_faf_threshold=bs1_supporting_faf_threshold,
        rare_variant_faf_threshold=rare_variant_faf_threshold,
        small_indel_size_threshold=small_indel_size_threshold,
        allele_count_threshold=allele_count_rare_variant_threshold,
        use_lcr=not no_lcr,
        missing_faf_suggests_absence=missing_faf_suggests_absence,
    )

    coverage_datasets = []
    coverage_meta = []   # parallel list of (gnomad_version, gnomad_data_type)
    for label, version, data_type, path in [
        ('gnomAD v4.1 joint', 'v4.1', 'joint',  coverage_v4_joint),
        ('gnomAD v4.1 exome', 'v4.1', 'exome',  coverage_v4_exome),
        ('gnomAD v3.1 genome', 'v3.1', 'genome', coverage_v3_genome),
    ]:
        if path is not None:
            print(f'Loading coverage: {path} ({label}) ...')
            coverage_datasets.append(read_coverage(path))
            coverage_meta.append((version, data_type))

    lcr_data = read_lcr(lcr) if lcr else None

    conn = psycopg2.connect(db_url, options=f'-c search_path={schema}')
    try:
        with conn.cursor() as cur:
            if vrs_digest:
                cur.execute(_QUERY_ONE, (vrs_digest,))
            elif overwrite:
                cur.execute(_QUERY_ALL)
            else:
                cur.execute(_QUERY_UNSCORED, (method_name,))
            rows = cur.fetchall()

        if vrs_digest and not rows:
            print(f'Error: VRS digest {vrs_digest!r} not found in the database.')
            return

        # Load every release's frequencies up front, so each variant can be
        # scored against the release whose coverage is sufficient for it.
        with conn.cursor() as cur:
            if vrs_digest:
                cur.execute(_QUERY_FREQ_ONE, (vrs_digest,))
            else:
                cur.execute(_QUERY_FREQ)
            freq_by_variant = {}
            for digest, version, data_type, flags, ac, faf, faf_pop in cur:
                freq_by_variant.setdefault(digest, {})[(version, data_type)] = (
                    flags, ac, faf, faf_pop)
        print(f'Loaded gnomAD frequencies for {len(freq_by_variant)} variant(s).')

        print(f'Computing evidence codes for {len(rows)} variant(s) ...')
        results = []
        for (vrs_digest_row, hgvs_cdna, chr_, pos, ref, alt) in rows:

            if debug:
                print(f'Analyzing {hgvs_cdna}')

            code, msg, dataset_idx = _compute_evidence_code(
                hgvs_cdna, chr_, pos, ref, alt,
                freq_by_variant.get(vrs_digest_row, {}),
                coverage_datasets, coverage_meta, lcr_data, config=config, debug=debug,
            )
            if dataset_idx is not None:
                gnomad_version, gnomad_data_type = coverage_meta[dataset_idx]
            else:
                gnomad_version, gnomad_data_type = None, None
            results.append((vrs_digest_row, code, msg, method_name, gnomad_version, gnomad_data_type))
            if dry_run:
                print(f'  {hgvs_cdna}  →  {code}')
                print(f'  {msg}')

        if dry_run:
            print(f'Dry run: {len(results)} result(s) computed, no database write.')
            return

        print(f'Writing {len(results)} rows to analysis_provisional_evidence_codes ...')
        with conn.cursor() as cur:
            for i in range(0, len(results), DB_BATCH):
                psycopg2.extras.execute_values(
                    cur,
                    """INSERT INTO analysis_provisional_evidence_codes
                           ("VRS_Digest", popfreq_code, popfreq_description, method_name,
                            gnomad_version, gnomad_data_type)
                       VALUES %s
                       ON CONFLICT ("VRS_Digest", method_name) DO UPDATE
                           SET popfreq_code        = EXCLUDED.popfreq_code,
                               popfreq_description = EXCLUDED.popfreq_description,
                               gnomad_version      = EXCLUDED.gnomad_version,
                               gnomad_data_type    = EXCLUDED.gnomad_data_type""",
                    results[i:i + DB_BATCH],
                )
        conn.commit()
        print(f'Done. Wrote {len(results)} rows.')
    finally:
        conn.close()


if __name__ == '__main__':
    main()
