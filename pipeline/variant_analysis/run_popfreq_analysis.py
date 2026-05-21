#!/usr/bin/env python
# coding: utf-8

"""
Populate analysis_provisional_evidence_codes from the database.

Reads gnomAD v4 data from report_gnomad, GRCh38 coordinates from
genomic_coordinates, and computes provisional population frequency
evidence codes for each variant. Results are upserted into
analysis_provisional_evidence_codes.

A gnomAD v4.1 coverage Parquet file and optionally a low-complexity
region BED file must still be supplied as local files.
"""

import bisect
import logging
import math
import os

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

SMALL_INDEL_SIZE_THRESHOLD = 50
ALLELE_COUNT_RARE_VARIANT_THRESHOLD = 1

BA1_FAF_THRESHOLD = 0.001
BS1_FAF_THRESHOLD = 0.0001
BS1_SUPPORTING_FAF_THRESHOLD = 0.00001
RARE_VARIANT_FAF_THRESHOLD = 0.00001

BA1_MSG = (f"The Total GrpMax filtering allele frequency (the lower threshold of the 95%% CI) "
           f"in gnomAD v4.1 is %s in the %s genetic ancestry group (based on %s/%s alleles) "
           f"which is above the ENIGMA BRCA1/2 VCEP threshold (>{BA1_FAF_THRESHOLD}) for BA1 (BA1 met).")
BS1_MSG = (f"The Total GrpMax filtering allele frequency (the lower threshold of the 95%% CI) "
           f"in gnomAD v4.1 is %s in the %s genetic ancestry group (based on %s/%s alleles) "
           f"which is above the ENIGMA BRCA1/2 VCEP threshold (>{BS1_FAF_THRESHOLD}) for BS1, "
           f"and below the BA1 threshold (>{BA1_FAF_THRESHOLD}) (BS1 met).")
BS1_SUPPORTING_MSG = (f"The Total GrpMax filtering allele frequency (the lower threshold of the 95%% CI) "
                      f"in gnomAD v4.1 is %s in the %s genetic ancestry group (based on %s/%s alleles) "
                      f"which is above the ENIGMA BRCA1/2 VCEP threshold (>{RARE_VARIANT_FAF_THRESHOLD}) "
                      f"for BS1_Supporting, and below the BS1 threshold (>{BS1_FAF_THRESHOLD}) (BS1_Supporting met).")
NO_CODE_MET_MSG = (f"The Total GrpMax filtering allele frequency (the lower threshold of the 95%% CI) "
                   f"in gnomAD v4.1 is %s in the %s genetic ancestry group (based on %s/%s alleles) "
                   f"which is below the ENIGMA BRCA1/2 VCEP threshold (>{BS1_SUPPORTING_FAF_THRESHOLD}) "
                   f"for BS1_Supporting and does not meet any population code "
                   f"(BA1, BS1, BS1_Supporting, PM2_Supporting are not met).")
NO_CODE_NO_FAF_MSG = (f"This variant is recorded in gnomAD v4.1, however the Total GrpMax filtering allele frequency "
                      f"(the lower threshold of the 95% CI) was not calculated, therefore this variant does not meet "
                      f"any population code (BA1, BS1, BS1_Supporting, PM2_Supporting are not met).")
NO_CODE_INDEL_MSG = (f"This [duplication/insertion/deletion/delins/large genomic rearrangement] variant "
                       f"was not observed in gnomAD v4.1, but PM2_Supporting was not applied since recall "
                       f"is considered suboptimal for this type of variant (PM2_Supporting not met). ")
PM2_SUPPORTING_MSG = "This variant is absent from gnomAD v4.1 (PM2_Supporting met)."
FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG_MSG = (f"This variant is present in gnomAD but is not meeting the "
                                                    f"specified read depths threshold ≥{READ_DEPTH_THRESHOLD_FREQUENT_VARIANT} "
                                                    f"(PM2_Supporting, BS1, and BA1 are not met).")
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


def analyze_one_dataset(faf95_popmax_str, allele_count, snv_or_small_indel,
                        read_depth, vcf_filter_flag,
                        allele_count_threshold,
                        chrom, genome_start, genome_end, lcr,
                        faf95_popmax, faf95_popmax_population,
                        allele_count_pop, allele_number_pop, debug=True):
    if vcf_filter_flag:
        return(FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG, FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG_MSG)
    rare_variant = False
    if field_defined(faf95_popmax_str):
        faf = float(faf95_popmax_str)
        if math.isnan(faf):
            rare_variant = True
        elif faf <= RARE_VARIANT_FAF_THRESHOLD:
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
        if faf > BA1_FAF_THRESHOLD:
            msg = BA1_MSG % (faf95_popmax, faf95_popmax_population, allele_count_pop, allele_number_pop)
            return(BA1, msg)
        elif faf > BS1_FAF_THRESHOLD:
            msg = BS1_MSG % (faf95_popmax, faf95_popmax_population, allele_count_pop, allele_number_pop)
            return(BS1, msg)
        elif faf > BS1_SUPPORTING_FAF_THRESHOLD:
            msg = BS1_SUPPORTING_MSG % (faf95_popmax, faf95_popmax_population, allele_count_pop, allele_number_pop)
            return(BS1_SUPPORTING, msg)
        else:
            msg = NO_CODE_MET_MSG % (faf95_popmax, faf95_popmax_population, allele_count_pop, allele_number_pop)
            return(NO_CODE, msg)
    if debug:
        print("Rare variant.  Allele count", allele_count, "SNV", snv_or_small_indel)
    if not field_defined(allele_count):
        meets_allele_count_threshold = False
    elif int(allele_count) <= allele_count_threshold:
        meets_allele_count_threshold = False
    else:
        meets_allele_count_threshold = True
    if meets_allele_count_threshold:
        if field_defined(faf95_popmax_str):
            return(NO_CODE, NO_CODE_NO_FAF_MSG)
        else:
            msg = NO_CODE_MET_MSG % (faf95_popmax, faf95_popmax_population, allele_count_pop, allele_number_pop)
            return(NO_CODE, msg)
    if lcr is not None:
        if overlaps_lcr(chrom, genome_start, genome_end, lcr):
            return(FAIL_LCR, FAIL_LCR_MSG)
    if snv_or_small_indel:
        return(PM2_SUPPORTING, PM2_SUPPORTING_MSG)
    else:
        return(NO_CODE_INDEL, NO_CODE_INDEL_MSG)


def _compute_evidence_code(hgvs_cdna, chr_, pos, ref,
                            gnomad_flags, gnomad_allele_count,
                            gnomad_faf95, gnomad_faf95_population,
                            cov4, lcr, debug=False):
    """Compute the popfreq evidence code for one variant from DB row fields."""
    start = int(pos)
    end = start + len(ref) - 1
    chrom = int(chr_.lstrip('chr'))

    _, read_depth = estimate_coverage(start, end, chrom, cov4, debug=debug)
    snv_or_small_indel = (end - start <= SMALL_INDEL_SIZE_THRESHOLD)

    # NULL gnomAD columns (variant absent from gnomAD v4) are treated as '-' (undefined).
    faf95 = gnomad_faf95 if gnomad_faf95 is not None else '-'
    allele_count = gnomad_allele_count if gnomad_allele_count is not None else '-'
    faf95_pop = gnomad_faf95_population if gnomad_faf95_population is not None else '-'
    # report_gnomad.Flags uses '-' for PASS; any other non-null value is a filter flag.
    is_flagged = gnomad_flags is not None and gnomad_flags != '-'

    if debug:
        print(f'  _compute_evidence_code inputs:')
        print(f'    hgvs_cdna={hgvs_cdna!r}  chr={chr_}  pos={pos}  ref={ref!r}')
        print(f'    gnomad_flags={gnomad_flags!r}  gnomad_allele_count={gnomad_allele_count!r}')
        print(f'    gnomad_faf95={gnomad_faf95!r}  gnomad_faf95_population={gnomad_faf95_population!r}')
        print(f'  derived:')
        print(f'    start={start}  end={end}  chrom={chrom}  snv_or_small_indel={snv_or_small_indel}')
        print(f'    faf95={faf95!r}  faf95_pop={faf95_pop!r}  allele_count={allele_count!r}')
        print(f'    is_flagged={is_flagged}  read_depth={read_depth}')

    # Per-population allele counts are not stored in the DB; '-' is used as placeholder
    # in the message text only (does not affect code assignment logic).
    code, msg = analyze_one_dataset(
        faf95, allele_count, snv_or_small_indel, read_depth, is_flagged,
        ALLELE_COUNT_RARE_VARIANT_THRESHOLD,
        chrom, start, end, lcr,
        faf95, faf95_pop,
        '-', '-',
        debug=debug,
    )

    if debug:
        print(f'  output: code={code!r}')

    return code, msg


_CREATE_TABLE = """
CREATE TABLE IF NOT EXISTS analysis_provisional_evidence_codes (
    "VRS_Digest"        text PRIMARY KEY
                             REFERENCES variant("VRS_Digest") ON DELETE CASCADE,
    popfreq_code        text,
    popfreq_description text
)
"""

_QUERY_ALL = """
SELECT v."VRS_Digest", v."HGVS_cDNA",
       gc.chr, gc.pos, gc.ref,
       rg."Flags", rg."Allele_count",
       rg.faf95_popmax, rg.faf95_popmax_population
FROM variant v
JOIN genomic_coordinates gc
     ON gc."VRS_Digest_id" = v."VRS_Digest" AND gc.assembly = 'GRCh38'
LEFT JOIN report_gnomad rg
     ON rg."VRS_Digest_id" = v."VRS_Digest" AND rg.version = 'v4'
"""

_QUERY_UNSCORED = _QUERY_ALL + """
LEFT JOIN analysis_provisional_evidence_codes apec
     ON apec."VRS_Digest" = v."VRS_Digest"
WHERE apec."VRS_Digest" IS NULL
"""

_QUERY_ONE = _QUERY_ALL + """
WHERE v."VRS_Digest" = %s
"""


@click.command()
@click.option('--db-url', default='postgresql://postgres:postgres@localhost/storage.pg',
              envvar='PIPELINE_DB_URL', show_default=True)
@click.option('--schema', default='pipeline', show_default=True)
@click.option('--coverage-file', required=True, type=click.Path(exists=True),
              help='Parquet file of gnomAD v4.1 per-position coverage (chrom, pos, mean)')
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
def main(db_url, schema, coverage_file, lcr, overwrite, debug, dry_run, vrs_digest):
    logging.basicConfig(level=logging.WARNING, format='%(levelname)s %(message)s')

    if debug:
        dry_run = True

    print(f'Loading coverage: {coverage_file} ...')
    cov4 = read_coverage(coverage_file)
    lcr_data = read_lcr(lcr) if lcr else None

    conn = psycopg2.connect(db_url, options=f'-c search_path={schema}')
    try:
        if not dry_run:
            with conn.cursor() as cur:
                cur.execute(_CREATE_TABLE)
            conn.commit()

        with conn.cursor() as cur:
            if vrs_digest:
                cur.execute(_QUERY_ONE, (vrs_digest,))
            elif overwrite:
                cur.execute(_QUERY_ALL)
            else:
                cur.execute(_QUERY_UNSCORED)
            rows = cur.fetchall()

        if vrs_digest and not rows:
            print(f'Error: VRS digest {vrs_digest!r} not found in the database.')
            return

        print(f'Computing evidence codes for {len(rows)} variant(s) ...')
        results = []
        for (vrs_digest_row, hgvs_cdna, chr_, pos, ref,
             flags, allele_count, faf95, faf95_pop) in rows:

            if debug:
                print(f'Analyzing {hgvs_cdna}')

            code, msg = _compute_evidence_code(
                hgvs_cdna, chr_, pos, ref,
                flags, allele_count, faf95, faf95_pop,
                cov4, lcr_data, debug=debug,
            )
            results.append((vrs_digest_row, code, msg))
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
                           ("VRS_Digest", popfreq_code, popfreq_description)
                       VALUES %s
                       ON CONFLICT ("VRS_Digest") DO UPDATE
                           SET popfreq_code        = EXCLUDED.popfreq_code,
                               popfreq_description = EXCLUDED.popfreq_description""",
                    results[i:i + DB_BATCH],
                )
        conn.commit()
        print(f'Done. Wrote {len(results)} rows.')
    finally:
        conn.close()


if __name__ == '__main__':
    main()
