#!/usr/bin/env python
# coding: utf-8

import bisect
import itertools
import argparse
import csv
from collections import OrderedDict
import math
import numpy as np
import pandas as pd
import sys

csv.field_size_limit(sys.maxsize)

POPFREQ_CODE_ID = "Provisional_Evidence_Code_Popfreq"
POPFREQ_CODE_DESCR = "Provisional_Evidence_Description_Popfreq"


BA1 = "BA1 (met)"
BS1 = "BS1 (met)"
BS1_SUPPORTING = "BS1_Supporting (met)"
NO_CODE = "No code met (below threshold)"
NO_CODE_NON_SNV = "No code met (non-SNV)"
PM2_SUPPORTING = "PM2_Supporting (met)"
FAIL_QC = "No code met (QC filter)"
FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG = "No code met (read depth, flags)"
FAIL_NEEDS_REVIEW = "No code met (needs review)"
FAIL_NEEDS_SOFTWARE_REVIEW = "No code met (needs software review)"
FAIL_LCR = "No code met (low-complexity region)"

READ_DEPTH_THRESHOLD_FREQUENT_VARIANT = 20
READ_DEPTH_THRESHOLD_RARE_VARIANT = 25

SMALL_INDEL_SIZE_THRESHOLD = 50
ALLELE_COUNT_RARE_VARIANT_THRESHOLD = 1
BA1_FAF_THRESHOLD = 0.001
BS1_FAF_THRESHOLD = 0.0001
BS1_SUPPORTING_FAF_THRESHOLD = 0.00001
RARE_VARIANT_FAF_THRESHOLD = 0.00001

# Messages used by analyze_across_datasets (v2/v3 wording, 2 format args: faf, population)
BA1_MSG = "The highest non-cancer, non-founder population filter allele frequency in gnomAD v2.1 (exomes only, non-cancer subset, read depth ≥20) or gnomAD v3.1 (non-cancer subset, read depth ≥20) is %s in the %s population, which is above the ENIGMA BRCA1/2 VCEP threshold (>0.001) for BA1 (BA1 met)."
BS1_MSG = "The highest non-cancer, non-founder population filter allele frequency in gnomAD v2.1 (exomes only, non-cancer subset, read depth ≥20) or gnomAD v3.1 (non-cancer subset, read depth ≥20) is %s in the %s population, which is above the ENIGMA BRCA1/2 VCEP threshold (>0.0001) for BS1, and below the BA1 threshold (>0.001) (BS1 met)."
BS1_SUPPORTING_MSG = "The highest non-cancer, non-founder population filter allele frequency in gnomAD v2.1 (exomes only, non-cancer subset, read depth ≥20) or gnomAD v3.1 f(non-cancer subset, read depth ≥20) is %s in the %s population which is within the ENIGMA BRCA1/2 VCEP threshold (>0.00002 to ≤ 0.0001) for BS1_Supporting (BS1_Supporting met)."
PM2_SUPPORTING_MSG = "This variant is absent from gnomAD v2.1 (exomes only, non-cancer subset, read depth ≥25) and gnomAD v3.1 (non-cancer subset, read depth ≥25) (PM2_Supporting met)."
NO_CODE_MSG = "This variant is present in gnomAD v2.1 (exomes only, non-cancer subset) or gnomAD v3.1 (non-cancer subset) but is below the ENIGMA BRCA1/2 VCEP threshold >0.00002 for BS1_Supporting (PM2_Supporting, BS1, and BA1 are not met)."
NO_CODE_NON_SNV_MSG = "This [insertion/deletion/large genomic rearrangement] variant was not observed in gnomAD v2.1 (exomes only, non-cancer subset) or gnomAD v3.1 (non-cancer subset), but PM2_Supporting was not applied since recall is suboptimal for this type of variant (PM2_Supporting not met)."
FAIL_QC_MSG = "This variant was not observed in gnomAD v2.1 (exomes only, non-cancer subset) or gnomAD v3.1 (non-cancer subset), but PM2_Supporting was not applied since this variant failed a gnomAD QC filter (PM2_Supporting not met)."
FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG_MSG = "This variant is present in gnomAD but is not meeting the specified read depths threshold ≥20 (PM2_Supporting, BS1, and BA1 are not met)."
FAIL_NEEDS_REVIEW_MSG = "No code is met (variant needs review)"
FAIL_NEEDS_SOFTWARE_REVIEW_MSG = "No code is met (variant needs software review)"
FAIL_LCR_MSG = "PM2_Supporting was not applied since this variant overlaps a low-complexity region (LCR)."

# Messages used by analyze_one_dataset (4 format args: faf, population, allele_count, allele_number)
BA1_MSG_V4 = (f"The Total GrpMax filtering allele frequency (the lower threshold of the 95%% CI) "
              f"in gnomAD v2.1/v3.1 is %s in the %s genetic ancestry group (based on %s/%s alleles) "
              f"which is above the ENIGMA BRCA1/2 VCEP threshold (>{BA1_FAF_THRESHOLD}) for BA1 (BA1 met).")
BS1_MSG_V4 = (f"The Total GrpMax filtering allele frequency (the lower threshold of the 95%% CI) "
              f"in gnomAD v2.1/v3.1 is %s in the %s genetic ancestry group (based on %s/%s alleles) "
              f"which is above the ENIGMA BRCA1/2 VCEP threshold (>{BS1_FAF_THRESHOLD}) for BS1, "
              f"and below the BA1 threshold (>{BA1_FAF_THRESHOLD}) (BS1 met).")
BS1_SUPPORTING_MSG_V4 = (f"The Total GrpMax filtering allele frequency (the lower threshold of the 95%% CI) "
                         f"in gnomAD v2.1/v3.1 is %s in the %s genetic ancestry group (based on %s/%s alleles) "
                         f"which is above the ENIGMA BRCA1/2 VCEP threshold (>{RARE_VARIANT_FAF_THRESHOLD}) "
                         f"for BS1_Supporting, and below the BS1 threshold (>{BS1_FAF_THRESHOLD}) (BS1_Supporting met).")
NO_CODE_MET_MSG = (f"The Total GrpMax filtering allele frequency (the lower threshold of the 95%% CI) "
                   f"in gnomAD v2.1/v3.1 is %s in the %s genetic ancestry group (based on %s/%s alleles) "
                   f"which is below the ENIGMA BRCA1/2 VCEP threshold (>{BS1_SUPPORTING_FAF_THRESHOLD}) "
                   f"for BS1_Supporting and does not meet any population code "
                   f"(BA1, BS1, BS1_Supporting, PM2_Supporting are not met).")
NO_CODE_NO_FAF_MSG = (f"This variant is recorded in gnomAD v2.1/v3.1, however the Total GrpMax filtering allele "
                      f"frequency (the lower threshold of the 95% CI) was not calculated, therefore this variant "
                      f"does not meet any population code (BA1, BS1, BS1_Supporting, PM2_Supporting are not met).")

MESSAGES_PER_CODE = {
    BA1: BA1_MSG,
    BS1: BS1_MSG,
    BS1_SUPPORTING: BS1_SUPPORTING_MSG,
    NO_CODE: NO_CODE_MSG,
    NO_CODE_NON_SNV: NO_CODE_NON_SNV_MSG,
    PM2_SUPPORTING: PM2_SUPPORTING_MSG,
    FAIL_QC: FAIL_QC_MSG,
    FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG: FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG_MSG,
    FAIL_NEEDS_REVIEW: FAIL_NEEDS_REVIEW_MSG,
    FAIL_NEEDS_SOFTWARE_REVIEW: FAIL_NEEDS_SOFTWARE_REVIEW_MSG,
    FAIL_LCR: FAIL_LCR_MSG
    }

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="built_final.tsv",
                        help="Input file with variant data")
    parser.add_argument("-o", "--output", default="built_with_popfreq.tsv",
                        help="Output file with new columns")
    parser.add_argument("-d", "--data_dir", default="./processed_brca",
                        help="Directory with the processed files")
    parser.add_argument("--popfreq-code-id", default=POPFREQ_CODE_ID,
                        help="Column header for the popfreq evidence code")
    parser.add_argument("--popfreq-code-descr", default=POPFREQ_CODE_DESCR,
                        help="Column header for the popfreq evidence description")
    parser.add_argument("--lcr", default=None,
                        help="BED file of low-complexity regions; variants overlapping an LCR "
                             "will not be assigned PM2_Supporting")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="Print debugging info")
    args = parser.parse_args()
    return(args)


def read_flags(flag_data):
    flags = {}
    for row in flag_data:
        flags[row["ID"]] = row
    return(flags)


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
    # Convert variant to 0-based half-open for comparison with BED
    var_start = start - 1
    var_end = end
    # Find the first region whose start >= var_end; none of those can overlap.
    # Search among region start coordinates using bisect.
    starts = [r[0] for r in regions]
    idx = bisect.bisect_left(starts, var_end)
    # Walk backwards from idx-1; stop as soon as a region ends before var_start.
    for i in range(idx - 1, -1, -1):
        r_start, r_end = regions[i]
        if r_end <= var_start:
            break
        return True  # r_start < var_end and r_end > var_start
    return False


def estimate_coverage(start, end, chrom, df_cov, debug=False, use_median=False):
    positions = list(range(start, end+1))
    coverage_this_chrom = df_cov.loc[df_cov["chrom"] == int(chrom)]
    positions_this_variant = coverage_this_chrom[coverage_this_chrom["pos"].isin(positions)]
    meanval = positions_this_variant["mean"].mean()
    medianval = positions_this_variant["median"].median()
    if pd.isna(meanval) or pd.isna(medianval):
        observable = False
        coverage = 0
    else:
        observable = True
        if use_median:
            coverage = min(meanval, medianval)
        else:
            coverage = meanval
    if debug:
        print("coverage assessment: observable:", observable, "meanval:",
              meanval, "medianval:", medianval, "coverage", coverage)
    return(observable, coverage)




def initialize_output_file(input_file, output_filename,
                           popfreq_code_id=POPFREQ_CODE_ID,
                           popfreq_code_descr=POPFREQ_CODE_DESCR):
    """
    Create an empty output file with the new columns
    """
    new_columns = [popfreq_code_id, popfreq_code_descr]
    input_header_row = input_file.fieldnames
    if "change_type" in input_header_row:
        idx = input_header_row.index("change_type")
        output_header_row = input_header_row[:idx] + new_columns \
            + input_header_row[idx:]
    else:
        output_header_row = input_header_row + new_columns
    output_file = csv.DictWriter(open(output_filename,"w"),
                                 fieldnames=output_header_row,
                                 delimiter = '\t')
    output_file.writeheader()
    return(output_file)


def field_defined(field):
    """
    Return a binary value indicating whether or not this variant has the popmax FAF defined
    """
    return(field != "-")


def get_population_allele_counts(population, variant, count_template, number_template):
    """
    Given a population code (e.g. "AMR"), return the population-specific allele count
    and allele number using the provided column name templates (with {} as placeholder).
    Returns (None, None) if population is None or "NA", or if the expected keys are absent.
    """
    if population is None or population == "NA":
        return None, None
    s1 = count_template.format(population)
    s2 = number_template.format(population)
    if s1 in variant and s2 in variant:
        return variant[s1], variant[s2]
    return None, None


def analyze_one_dataset(faf95_popmax_str, allele_count, snv_or_small_indel,
                        read_depth, vcf_filter_flag,
                        allele_count_threshold,
                        chrom, genome_start, genome_end, lcr,
                        faf95_popmax, faf95_popmax_population,
                        allele_count_pop, allele_number_pop, debug=True):
    #
    # Rule out error conditions: VCF filter flag, low coverage.
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
    #
    # Address the cases where FAF is defined, and the variant is a candidate for a
    # evidence code for high population frequency (BA1, BS1, BS1_SUPPORTING)
    if not rare_variant:
        if debug:
            print("Not rare variant.  FAF:", faf)
        if faf > BA1_FAF_THRESHOLD:
            msg = BA1_MSG_V4 % (faf95_popmax, faf95_popmax_population, allele_count_pop, allele_number_pop)
            return(BA1, msg)
        elif faf > BS1_FAF_THRESHOLD:
            msg = BS1_MSG_V4 % (faf95_popmax, faf95_popmax_population, allele_count_pop, allele_number_pop)
            return(BS1, msg)
        elif faf > BS1_SUPPORTING_FAF_THRESHOLD:
            msg = BS1_SUPPORTING_MSG_V4 % (faf95_popmax, faf95_popmax_population, allele_count_pop, allele_number_pop)
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
    #
    # If it meets the allele count threshold, then it's a no-code-met variant.  Choose the message by
    # whether or not the FAF has been estimated.
    if meets_allele_count_threshold:
        if field_defined(faf95_popmax_str):
            return(NO_CODE, NO_CODE_NO_FAF_MSG)
        else:
            msg = NO_CODE_MET_MSG % (faf95_popmax, faf95_popmax_population, allele_count_pop, allele_number_pop)
            return(NO_CODE, msg)
    #
    # At this point, it's PM2_Supporting, unless it overlaps with an LCR or it is a non-SNV.
    if lcr is not None:
        if overlaps_lcr(chrom, genome_start, genome_end, lcr):
            return(FAIL_LCR, FAIL_LCR_MSG)
    if snv_or_small_indel:
        return(PM2_SUPPORTING, PM2_SUPPORTING_MSG)
    else:
        return(NO_CODE_NON_SNV, NO_CODE_NON_SNV_MSG)


def analyze_across_datasets(code_v2, faf_v2, faf_popmax_v2, in_v2,
                            code_v3, faf_v3, faf_popmax_v3, in_v3, is_snv,
                            debug=False):
    """
    Given the per-dataset evidence codes, generate an overall evidence code
    """
    benign_codes = [BA1, BS1, BS1_SUPPORTING]
    pathogenic_codes = [PM2_SUPPORTING]
    intermediate_codes = [ NO_CODE, NO_CODE_NON_SNV]
    fail_codes = [FAIL_QC, FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG,
                  FAIL_NEEDS_REVIEW, FAIL_NEEDS_SOFTWARE_REVIEW,
                  FAIL_LCR]
    ordered_success_codes = benign_codes + intermediate_codes + pathogenic_codes
    
    success_codes = set(ordered_success_codes)
    failure_codes = set(fail_codes)
    ordered_codes = ordered_success_codes + list(failure_codes)

    #
    # First, rule out the case of outright contradictions
    if code_v2 in benign_codes and code_v3 in pathogenic_codes \
       or code_v2 in pathogenic_codes and code_v3 in benign_codes:
        return(FAIL_NEEDS_REVIEW, FAIL_NEEDS_REVIEW_MSG)
    #
    # Next, rule out the case where neither dataset has reliable data.
    # Arbitrarily use the message for v3, as the newer and presumably more robust.
    if code_v2 in failure_codes and code_v3 in failure_codes:
        return(code_v3, MESSAGES_PER_CODE[code_v3])

    #
    # At this point, we can assume that at least one dataset has
    # reliable data.
    #
    # Next, if both datasets point to a pathogenic effect, or
    # if one points to a pathogenic effect and the other an error,
    # then return the pathogenic effect.
    if code_v2 in pathogenic_codes and code_v3 in pathogenic_codes:
        return(code_v2, MESSAGES_PER_CODE[code_v2])
    elif code_v2 in pathogenic_codes and code_v3 in failure_codes:
        return(code_v2, MESSAGES_PER_CODE[code_v2])
    elif code_v3 in pathogenic_codes and code_v2 in failure_codes:
        return(code_v3, MESSAGES_PER_CODE[code_v3])

    #
    # Next, if either dataset has an intermediate effect and the other
    # is not stronger (i.e. is also intermediate, or pathogenic, or failure),
    # return the intermediate effect code.
    if code_v2 in intermediate_codes and ordered_codes.index(code_v2) <= ordered_codes.index(code_v3):
        return(code_v2, MESSAGES_PER_CODE[code_v2])
    elif code_v3 in intermediate_codes and ordered_codes.index(code_v3) <= ordered_codes.index(code_v2):
        return(code_v3, MESSAGES_PER_CODE[code_v3])

    #
    # Now, at least one dataset must have a success code.  We can also assume that
    # neither is pathogenic (i.e. boht are benign, intermediate or failure).
    # In this case, identify and return the stronger code.
    if debug:
        print("prior to assertions, codes are", code_v2, code_v3)
    assert(code_v2 in benign_codes or code_v2 in intermediate_codes or code_v3 in benign_codes or code_v3 in intermediate_codes)
    if code_v3 == BA1:
        return(BA1, BA1_MSG % (faf_v3, faf_popmax_v3))
    elif code_v2 == BA1:
        return(BA1, BA1_MSG % (faf_v2, faf_popmax_v2))
    elif code_v3 == BS1:
        return(BS1, BS1_MSG % (faf_v3, faf_popmax_v3))
    elif code_v2 == BS1:
        return(BS1, BS1_MSG % (faf_v2, faf_popmax_v2))
    elif code_v3 == BS1_SUPPORTING:
        return(BS1_SUPPORTING, BS1_SUPPORTING_MSG % (faf_v3, faf_popmax_v3))
    elif code_v2 == BS1_SUPPORTING:
        return(BS1_SUPPORTING, BS1_SUPPORTING_MSG % (faf_v2, faf_popmax_v2))
    elif code_v2 == NO_CODE:
        return(code_v2, NO_CODE_MSG)
    elif code_v3 == NO_CODE:
        return(code_v3, NO_CODE_MSG)
    #
    # If we get here, there is a hole in the logic
    return(FAIL_NEEDS_REVIEW, FAIL_NEEDS_SOFTWARE_REVIEW_MSG)


def variant_is_flagged(variant_id, flags):
    if not field_defined(variant_id):
        return False
    assert(variant_id in flags)
    variant_flagged = False
    if flags[variant_id]["Filters"] != "PASS":
        variant_flagged = True
    return(variant_flagged)


def analyze_variant(variant, coverage_v2, coverage_v3, flags_v2, flags_v3,
                    popfreq_code_id=POPFREQ_CODE_ID,
                    popfreq_code_descr=POPFREQ_CODE_DESCR,
                    lcr=None, debug=False):
    """
    Analyze a single variant, adding the output columns
    """
    (observable_in_v2,
     read_depth_v2) = estimate_coverage(int(variant["pyhgvs_Hg37_Start"]),
                                        int(variant["pyhgvs_Hg37_End"]),
                                        int(variant["Chr"]),coverage_v2,
                                        debug=debug)
    (observable_in_v3,
     read_depth_v3) = estimate_coverage(int(variant["Hg38_Start"]),
                                        int(variant["Hg38_End"]),
                                        int(variant["Chr"]),coverage_v3,
                                        debug=debug)
    is_snv = (variant["Hg38_Start"] == variant["Hg38_End"]
              and len(variant["Ref"]) == 1 and len(variant["Alt"]) == 1)
    snv_or_small_indel = (is_snv or
                          (int(variant["Hg38_End"]) - int(variant["Hg38_Start"]) <= SMALL_INDEL_SIZE_THRESHOLD))
    if snv_or_small_indel:
        if observable_in_v2:
            variant_v2_code_id = PM2_SUPPORTING
        else:
            variant_v2_code_id = FAIL_QC
        if observable_in_v3:
            variant_v3_code_id = PM2_SUPPORTING
        else:
            variant_v3_code_id = FAIL_QC
    else:
        variant_v2_code_id = NO_CODE_NON_SNV
        variant_v3_code_id = NO_CODE_NON_SNV

    if debug:
        print("variant", variant["pyhgvs_cDNA"], " is snv:", is_snv, "preliminary codes", variant_v2_code_id,
              variant_v3_code_id, "observable", observable_in_v2, observable_in_v3)

    variant_in_v2 = False
    if (field_defined(variant["Variant_id_GnomAD"])
        and field_defined(variant["Allele_count_exome_GnomAD"])
        and observable_in_v2):
        variant_in_v2 = True
        population_v2 = variant["faf95_popmax_population_exome_GnomAD"]
        allele_count_pop_v2, allele_number_pop_v2 = get_population_allele_counts(
            population_v2, variant,
            "Allele_count_exome_{}_GnomAD", "Allele_number_exome_{}_GnomAD")
        (variant_v2_code_id, _) = analyze_one_dataset(
            variant["faf95_popmax_exome_GnomAD"],
            variant["Allele_count_exome_GnomAD"],
            snv_or_small_indel, read_depth_v2,
            variant_is_flagged(variant["Variant_id_GnomAD"], flags_v2),
            ALLELE_COUNT_RARE_VARIANT_THRESHOLD,
            int(variant["Chr"]), int(variant["pyhgvs_Hg37_Start"]), int(variant["pyhgvs_Hg37_End"]),
            None,  # no LCR for hg19 coords
            variant["faf95_popmax_exome_GnomAD"],
            population_v2,
            allele_count_pop_v2, allele_number_pop_v2,
            debug=debug)
        if debug:
            print("gnomAD V2 variant", variant["Variant_id_GnomAD"],
                  "popmax", variant["faf95_popmax_exome_GnomAD"],
                  "allele count", variant["Allele_count_exome_GnomAD"],
                  "read depth", read_depth_v2, "snv", is_snv,
                  "V2 code", variant_v2_code_id)
    variant_in_v3 = False
    if (field_defined(variant["Variant_id_GnomADv3"]) and observable_in_v3):
        variant_in_v3 = True
        population_v3 = variant["faf95_popmax_population_genome_GnomADv3"]
        allele_count_pop_v3, allele_number_pop_v3 = get_population_allele_counts(
            population_v3, variant,
            "Allele_count_genome_{}_GnomADv3", "Allele_number_genome_{}_GnomADv3")
        (variant_v3_code_id, _) = analyze_one_dataset(
            variant["faf95_popmax_genome_GnomADv3"],
            variant["Allele_count_genome_GnomADv3"],
            snv_or_small_indel, read_depth_v3,
            variant_is_flagged(variant["Variant_id_GnomADv3"], flags_v3),
            ALLELE_COUNT_RARE_VARIANT_THRESHOLD,
            int(variant["Chr"]), int(variant["Hg38_Start"]), int(variant["Hg38_End"]),
            lcr,
            variant["faf95_popmax_genome_GnomADv3"],
            population_v3,
            allele_count_pop_v3, allele_number_pop_v3,
            debug=debug)
        if debug:
            print("gnomAD V3 variant", variant["Variant_id_GnomADv3"],
                  "popmax", variant["faf95_popmax_genome_GnomADv3"],
                  "allele count", variant["Allele_count_genome_GnomADv3"],
                  "read depth", read_depth_v3, "snv", is_snv,
                  "V3 code", variant_v3_code_id)
    (variant[popfreq_code_id],
     variant[popfreq_code_descr]) = analyze_across_datasets(variant_v2_code_id,variant["faf95_popmax_exome_GnomAD"],
                                                            variant["faf95_popmax_population_exome_GnomAD"],
                                                            variant_in_v2, variant_v3_code_id,
                                                            variant["faf95_popmax_genome_GnomADv3"],
                                                            variant["faf95_popmax_population_genome_GnomADv3"],
                                                            variant_in_v3, is_snv, debug)
    if debug:
        print("variant", variant["pyhgvs_cDNA"], "consensus code:", variant[popfreq_code_id], "msg",
              variant[popfreq_code_descr], "v2 code", variant_v2_code_id,
              "v3 code", variant_v3_code_id)
    return()



def main():
    args = parse_args()
    df_cov2 = pd.read_parquet(args.data_dir + '/df_cov_v2.parquet')
    df_cov3 = pd.read_parquet(args.data_dir + '/df_cov_v3.parquet')
    flags_v2 = read_flags(csv.DictReader(open(args.data_dir + "/brca.gnomAD.2.1.1.hg19.flags.tsv"),
                                         delimiter = "\t"))
    flags_v3 = read_flags(csv.DictReader(open(args.data_dir + "/brca.gnomAD.3.1.1.hg38.flags.tsv"),
                                         delimiter = "\t"))
    lcr = read_lcr(args.lcr) if args.lcr else None
    input_file = csv.DictReader(open(args.input), delimiter = "\t")
    output_file = initialize_output_file(input_file, args.output,
                                         args.popfreq_code_id,
                                         args.popfreq_code_descr)
    for variant in input_file:
        if args.debug:
            print("About to analyze", variant["pyhgvs_cDNA"])
        analyze_variant(variant, df_cov2, df_cov3, flags_v2, flags_v3,
                        popfreq_code_id=args.popfreq_code_id,
                        popfreq_code_descr=args.popfreq_code_descr,
                        lcr=lcr, debug=args.debug)
        output_file.writerow(variant)


if __name__ == "__main__":
    main()
