#!/usr/bin/env python
# coding: utf-8

import itertools
import argparse
import csv
from collections import OrderedDict
import math
import numpy as np
import pandas as pd

#import gnomad.variant_scoring.constants as cnts
#from common import config, hgvs_utils, variant_utils, utils
#from data_merging.brca_pseudonym_generator import _normalize_genomic_fnc

POPFREQ_CODE_ID = "Provisional_Evidence_Code_Popfreq"
POPFREQ_CODE_DESCR = "Provisional_Evidence_Description_Popfreq"


BA1 = "BA1 (met)"
BS1 = "BS1 (met)"
BS1_SUPPORTING = "BS1_Supporting (met)"
NO_CODE = "No code met (below threshold)"
NO_CODE_NON_SNV = "No code met for population data (indel)"
PM2_SUPPORTING = "PM2_Supporting (met)"
FAIL_NOT_ASSAYED = "Fail_Not_Assayed"
FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG = "No code met (read depth, flags)"
FAIL_NEEDS_REVIEW = "No code met (needs review)"
FAIL_NEEDS_SOFTWARE_REVIEW = "No code met (needs software review)"

READ_DEPTH_THRESHOLD_FREQUENT_VARIANT = 20
READ_DEPTH_THRESHOLD_RARE_VARIANT = 25

BA1_MSG = "The highest non-cancer, non-founder population filter allele frequency in gnomAD v2.1 (exomes only, non-cancer subset, read depth ≥20) or gnomAD v3.1 (non-cancer subset, read depth ≥20) is %s in the %s population, which is above the ENIGMA BRCA1/2 VCEP threshold (>0.001) for BA1 (BA1 met)."
BS1_MSG = "The highest non-cancer, non-founder population filter allele frequency in gnomAD v2.1 (exomes only, non-cancer subset, read depth ≥20) or gnomAD v3.1 (non-cancer subset, read depth ≥20) is %s in the %s population, which is above the ENIGMA BRCA1/2 VCEP threshold (>0.0001) for BS1, and below the BA1 threshold (>0.001) (BS1 met)."
BS1_SUPPORTING_MSG = "The highest non-cancer, non-founder population filter allele frequency in gnomAD v2.1 (exomes only, non-cancer subset, read depth ≥20) or gnomAD v3.1 f(non-cancer subset, read depth ≥20) is %s in the %s population which is within the ENIGMA BRCA1/2 VCEP threshold (>0.00002 to ≤ 0.0001) for BS1_Supporting (BS1_Supporting met)."
PM2_SUPPORTING_MSG = "This variant is absent from gnomAD v2.1 (exomes only, non-cancer subset, read depth ≥25) and gnomAD v3.1 (non-cancer subset, read depth ≥25) (PM2_Supporting met)."
NO_CODE_MSG = "This variant is present in gnomAD v2.1 (exomes only, non-cancer subset) or gnomAD v3.1 (non-cancer subset) but is below the ENIGMA BRCA1/2 VCEP threshold >0.00002 for BS1_Supporting (PM2_Supporting, BS1, and BA1 are not met)."
NO_CODE_NON_SNV_MSG = "This [insertion/deletion/large genomic rearrangement] variant was not observed in gnomAD v2.1 (exomes only, non-cancer subset) or gnomAD v3.1 (non-cancer subset), but PM2_Supporting was not applied since recall is suboptimal for this type of variant (PM2_Supporting not met)."
FAIL_NOT_ASSAYED_MSG = "Variant not tested in this dataset"
FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG_MSG = "This variant is present in gnomAD v2.1 (exomes only, non-cancer subset) or gnomAD v3.1 (non-cancer subset) but is not meeting the specified read depths threshold ≥20 OR was flagged as suspect by gnomAD (PM2_Supporting, BS1, and BA1 are not met)."
FAIL_NEEDS_REVIEW_MSG = "No code is met (variant needs review)"
FAIL_NEEDS_SOFTWARE_REVIEW_MSG = "No code is met (variant needs software review)"

MESSAGES_PER_CODE = {
    BA1: BA1_MSG,
    BS1: BS1_MSG,
    BS1_SUPPORTING: BS1_SUPPORTING_MSG,
    NO_CODE: NO_CODE_MSG,
    NO_CODE_NON_SNV: NO_CODE_NON_SNV_MSG,
    PM2_SUPPORTING: PM2_SUPPORTING_MSG,
    FAIL_NOT_ASSAYED: FAIL_NOT_ASSAYED_MSG,
    FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG: FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG_MSG,
    FAIL_NEEDS_REVIEW: FAIL_NEEDS_REVIEW_MSG,
    FAIL_NEEDS_SOFTWARE_REVIEW: FAIL_NEEDS_SOFTWARE_REVIEW_MSG
    }

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="built_final.tsv",
                        help="Input file with variant data")
    parser.add_argument("-o", "--output", default="built_with_popfreq.tsv",
                        help="Output file with new columns")
    parser.add_argument("-d", "--data_dir", default="./processed_brca",
                        help="Directory with the processed files")
    args = parser.parse_args()
    return(args)


def read_flags(flag_data):
    flags = {}
    for row in flag_data:
        flags[row["ID"]] = row
    return(flags)
    

def estimate_coverage(start, end, chrom, df_cov):
    positions = list(range(start, end+1))
    coverage_this_chrom = df_cov.loc[df_cov["chrom"] == int(chrom)]
    positions_this_variant = coverage_this_chrom[coverage_this_chrom["pos"].isin(positions)]
    meanval = positions_this_variant["mean"].mean()
    medianval = positions_this_variant["median"].median()
    return(meanval, medianval)




def initialize_output_file(input_file, output_filename):
    """
    Create an empty output file with the new columns
    """
    new_columns = [POPFREQ_CODE_ID, POPFREQ_CODE_DESCR]
    input_header_row = input_file.fieldnames
    if input_header_row.index("change_type"):
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


def is_variant_observable(start, end, chrom, coverage):
    #
    # Read the coverage at the variant position.  If the coverage metrics returned
    # are not a number, that means that the variant could not be tested in the
    # indicated dataset.
    (mean_read_depth, median_read_depth) = estimate_coverage(int(start),int(end),
                                                             int(chrom),coverage)
    if pd.isna(mean_read_depth) and pd.isna(median_read_depth):
        return(False)
    else:
        return(True)
    

def analyze_one_dataset(faf95_popmax_str, faf95_population, allele_count, is_snv,
                        read_depth, vcf_filter_flag, debug=True):
    #
    # Get the coverage data.  Rule out error conditions: low coverage, VCF filter flag.
    rare_variant = False
    if field_defined(faf95_popmax_str):
        faf = float(faf95_popmax_str)
        if pd.isna(faf):
            rare_variant = True
        elif faf <= 0.00002:
            rare_variant = True
    else:
        faf = None
        rare_variant = True
    if rare_variant and read_depth < READ_DEPTH_THRESHOLD_RARE_VARIANT:
        return(FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG)
    if (not rare_variant) and read_depth < READ_DEPTH_THRESHOLD_FREQUENT_VARIANT:
        return(FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG)
    if vcf_filter_flag:
        return(FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG)
    #
    # Address the cases where FAF is defined, and the variant is a candidate for a
    # evidence code for high population frequency (BA1, BS1, BS1_SUPPORTING)
    if not rare_variant:
        if faf > 0.001:
            return(BA1)
        elif faf >  0.0001:
            return(BS1)
        elif faf > 0.00002:
            return(BS1_SUPPORTING)
    if rare_variant:
        if not field_defined(allele_count):
            return(NO_CODE)
        elif int(allele_count) > 0:
            return(NO_CODE)
        else:
            if is_snv:
                return(PM2_SUPPORTING)
            else:
                return(NO_CODE_NON_SNV)
    return(NEEDS_REVIEW)


def analyze_across_datasets(code_v2, faf_v2, faf_popmax_v2, in_v2,
                            code_v3, faf_v3, faf_popmax_v3, in_v3,
                            debug=True):
    """
    Given the per-dataset evidence codes, generate an overall evidence code
    """
    benign_codes = [BA1, BS1, BS1_SUPPORTING]
    pathogenic_codes = [PM2_SUPPORTING]
    intermediate_codes = [ NO_CODE, NO_CODE_NON_SNV]
    ordered_success_codes = benign_codes + intermediate_codes + pathogenic_codes
    success_codes = set(ordered_success_codes)
    failure_codes = set([FAIL_NOT_ASSAYED, FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG])
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
    assert(variant_id in flags)
    variant_flagged = False
    if flags[variant_id]["Filters"] != "PASS":
        variant_flagged = True
    return(variant_flagged)


def analyze_variant(variant, coverage_v2, coverage_v3, flags_v2, flags_v3,
                    debug=True):
    """
    Analyze a single variant, adding the output columns
    """
    # Initialize the output columns.  First, check the read depth.  If the
    # read depth is defined at the variant position, this means that the
    # variant could have been observed in principle; set the default to
    # PM2_SUPPORTING (indicating that by default, the variant could have been
    # observed, but wasn't).  If no read depth is defined, this means that
    # the variant could not have been observed in the dataset in question
    # (for example, a deep intronic variant in an exome dataset).
    # In such a case, set the default value to NOT_ASSAYED.
    variant_v2_code_id = FAIL_NOT_ASSAYED
    variant_v3_code_id = FAIL_NOT_ASSAYED
    variant[POPFREQ_CODE_ID] = FAIL_NOT_ASSAYED
    variant[POPFREQ_CODE_DESCR] = FAIL_NOT_ASSAYED_MSG
    observable_in_v2 = False
    observable_in_v3 = False
    if is_variant_observable(int(variant["pyhgvs_Hg37_Start"]),int(variant["pyhgvs_Hg37_End"]),
                             int(variant["Chr"]),coverage_v2):
        variant_v2_code_id = PM2_SUPPORTING
        variant[POPFREQ_CODE_ID] = PM2_SUPPORTING
        variant[POPFREQ_CODE_DESCR] = PM2_SUPPORTING_MSG
        observable_in_v2 = True
    if is_variant_observable(int(variant["Hg38_Start"]),
                             int(variant["Hg38_End"]),
                             int(variant["Chr"]), coverage_v3):
        variant_v3_code_id = PM2_SUPPORTING
        variant[POPFREQ_CODE_ID] = PM2_SUPPORTING
        variant[POPFREQ_CODE_DESCR] = PM2_SUPPORTING_MSG
        observable_in_v3 = True
    is_snv = (variant["Hg38_Start"] == variant["Hg38_End"])
    if debug:
        print("variant is snv:", is_snv)
    #
    # the gnomAD v2 variant ID is set when the variant is in the genome
    # OR exome portion of gnomAD. Focus on variants that are in the exome
    # data by making sure that the allele count is defined.  The allele
    # count is the total number of observations of the variant in the gnomAD
    # dataset
    variant_in_v2 = False
    if (field_defined(variant["Variant_id_GnomAD"])
        and field_defined(variant["Allele_count_exome_GnomAD"])
        and observable_in_v2):
        variant_in_v2 = True 
        (mean_read_depth,
         median_read_depth) = estimate_coverage(int(variant["pyhgvs_Hg37_Start"]),
                                                int(variant["pyhgvs_Hg37_End"]),
                                                int(variant["Chr"]),coverage_v2)
        read_depth_v2 = min(mean_read_depth, median_read_depth)
        if pd.isna(read_depth_v2):
            read_depth_v2 = 0
        variant_v2_code_id = analyze_one_dataset(variant["faf95_popmax_exome_GnomAD"],
                                                         variant["faf95_popmax_population_exome_GnomAD"],
                                                         variant["Allele_count_exome_GnomAD"],
                                                         is_snv, read_depth_v2,
                                                         variant_is_flagged(variant["Variant_id_GnomAD"],
                                                                            flags_v2))
        if debug:
            print("gnomAD V2 variant", variant["Variant_id_GnomAD"],
                  "popmax", variant["faf95_popmax_exome_GnomAD"],
                  "allele count", variant["Allele_count_exome_GnomAD"],
                  "mean read depth", mean_read_depth,
                  "median read depth", median_read_depth, "snv", is_snv,
                  "V2 code", variant_v2_code_id)
    variant_in_v3 = False
    if (field_defined(variant["Variant_id_GnomADv3"]) and observable_in_v3):
        variant_in_v3 = True
        (mean_read_depth,
         median_read_depth) = estimate_coverage(int(variant["Hg38_Start"]),
                                                int(variant["Hg38_End"]),
                                                int(variant["Chr"]),
                                                coverage_v3)
        read_depth_v3 = min(mean_read_depth, median_read_depth)
        if pd.isna(read_depth_v3):
            read_depth_v3 = 0
        variant_v3_code_id = analyze_one_dataset(variant["faf95_popmax_genome_GnomADv3"],
                                                         variant["faf95_popmax_population_genome_GnomADv3"],
                                                         variant["Allele_count_genome_GnomADv3"],
                                                         is_snv, read_depth_v3,
                                                         variant_is_flagged(variant["Variant_id_GnomADv3"],
                                                                            flags_v3))
        if debug:
            print("gnomAD V3 variant", variant["Variant_id_GnomADv3"],
                  "popmax", variant["faf95_popmax_genome_GnomADv3"],
                  "allele count", variant["Allele_count_genome_GnomADv3"],
                  "mean read depth", mean_read_depth,
                  "median read depth", median_read_depth, "snv", is_snv,
                  "V3 code", variant_v3_code_id)
    (variant[POPFREQ_CODE_ID],
     variant[POPFREQ_CODE_DESCR]) = analyze_across_datasets(variant_v2_code_id,variant["faf95_popmax_exome_GnomAD"],
                                                            variant["faf95_popmax_population_exome_GnomAD"],
                                                            variant_in_v2, variant_v3_code_id, 
                                                            variant["faf95_popmax_genome_GnomADv3"],
                                                            variant["faf95_popmax_population_genome_GnomADv3"],
                                                            variant_in_v3)
    if debug:
        print("consensus code:", variant[POPFREQ_CODE_ID], "msg",
              variant[POPFREQ_CODE_DESCR])
    return()



def main():
    args = parse_args()
    df_cov2 = pd.read_parquet(args.data_dir + '/df_cov_v2.parquet')
    df_cov3 = pd.read_parquet(args.data_dir + '/df_cov_v3.parquet')
    flags_v2 = read_flags(csv.DictReader(open(args.data_dir + "/brca.gnomAD.2.1.1.hg19.flags.tsv"),
                                         delimiter = "\t"))
    flags_v3 = read_flags(csv.DictReader(open(args.data_dir + "/brca.gnomAD.3.1.1.hg38.flags.tsv"),
                                         delimiter = "\t"))
    input_file = csv.DictReader(open(args.input), delimiter = "\t")
    output_file = initialize_output_file(input_file, args.output)
    for variant in input_file:
        analyze_variant(variant, df_cov2, df_cov3, flags_v2, flags_v3)
        output_file.writerow(variant)


if __name__ == "__main__":
    main()

