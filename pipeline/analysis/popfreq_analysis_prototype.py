#!/usr/bin/env python
# coding: utf-8

import argparse
import csv
import math
import re
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

READ_DEPTH_THRESHOLD_FREQUENT_VARIANT = 20
READ_DEPTH_THRESHOLD_RARE_VARIANT = 25

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
    FAIL_NEEDS_SOFTWARE_REVIEW: FAIL_NEEDS_SOFTWARE_REVIEW_MSG
    }

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="test_build.tsv",
                        help="Input file with variant data")
    parser.add_argument("-o", "--output", default="test_out.tsv",
                        help="Output file with new columns")
    parser.add_argument("-d", "--data_dir", default="/Users/melissacline/Desktop/gnomAD/output",
                        help="Directory with the processed files")
#    parser.add_argument("-i", "--input", default="built_final.tsv",
#                        help="Input file with variant data")
#    parser.add_argument("-o", "--output", default="built_with_popfreq.tsv",
#                        help="Output file with new columns")
#    parser.add_argument("-d", "--data_dir", default="./processed_brca",
#                        help="Directory with the processed files")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="Print debugging info")
    args = parser.parse_args()
    return(args)


def read_flags(flag_data):
    flags = {}
    for row in flag_data:
        flags[row["ID"]] = row
    return(flags)


def read_coverage(coverage_file):
    """
    Read coverage data from CSV file and organize by chromosome and position.
    Returns: dict with structure {chrom: {pos: {"mean": float, "median": float}}}
    """
    coverage = {}
    reader = csv.DictReader(open(coverage_file), delimiter=',')
    for row in reader:
        chrom = int(row["chrom"])
        if chrom not in coverage:
            coverage[chrom] = {}
        pos = int(row["pos"])
        mean = float(row["mean"])
        if "median" in row:
            median = float(row["median"])
            coverage[chrom][pos] = {"mean": mean, "median": median}
        else:
            coverage[chrom][pos] = {"mean": mean}

    return coverage


def estimate_coverage(start, end, chrom, cov_data, debug=False, use_median=False):
    """
    Estimate coverage for a genomic region from coverage data stored as a dictionary.
    cov_data: dict with structure {chrom: {pos: {"mean": float, "median": float}}}
    """
    positions = list(range(start, end+1))
    chrom_key = int(chrom)

    # Get coverage values for positions in this variant
    mean_values = []
    median_values = []

    if chrom_key in cov_data:
        for pos in positions:
            if pos in cov_data[chrom_key]:
                mean_values.append(cov_data[chrom_key][pos]["mean"])
                if use_median:
                    median_values.append(cov_data[chrom_key][pos]["median"])

    # Calculate mean of means, and median of medians if applicable
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




def initialize_output_file(input_file, output_filename):
    """
    Create an empty output file with the new columns
    """
    new_columns = [POPFREQ_CODE_ID, POPFREQ_CODE_DESCR, "V4_Popfreq"]
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

    

def analyze_one_dataset(faf95_popmax_str, faf95_population, allele_count, is_snv,
                        read_depth, vcf_filter_flag, debug=True, test_coverage=True):
    #
    # Get the coverage data.  Rule out error conditions: low coverage, VCF filter flag.
    rare_variant = False
    if field_defined(faf95_popmax_str):
        faf = float(faf95_popmax_str)
        if math.isnan(faf):
            rare_variant = True
        elif faf <= 0.00002:
            rare_variant = True
    else:
        faf = None
        rare_variant = True
    if debug:
        print("Rare variant", rare_variant, "read depth", read_depth)
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
        if debug:
            print("Not rare variant.  FAF:", faf)
        if faf > 0.001:
            return(BA1)
        elif faf >  0.0001:
            return(BS1)
        elif faf > 0.00002:
            return(BS1_SUPPORTING)
    if rare_variant:
        if debug:
            print("Rare variant.  Allele count", allele_count, "SNV", is_snv)
        if not field_defined(allele_count):
            return(NO_CODE)
        elif int(allele_count) > 0:
            return(NO_CODE)
        else:
            if debug:
                print("Returning PM2_supporting or no code: is_snv", is_snv)
            if is_snv:
                return(PM2_SUPPORTING)
            else:
                return(NO_CODE_NON_SNV)
    return(NEEDS_REVIEW)


def analyze_across_datasets(code_v2, faf_v2, faf_popmax_v2, in_v2,
                            code_v3, faf_v3, faf_popmax_v3, in_v3, is_snv,
                            debug=False):
    """
    Given the per-dataset evidence codes, generate an overall evidence code
    """
    benign_codes = [BA1, BS1, BS1_SUPPORTING]
    pathogenic_codes = [PM2_SUPPORTING]
    intermediate_codes = [ NO_CODE, NO_CODE_NON_SNV]
    ordered_success_codes = benign_codes + intermediate_codes + pathogenic_codes
    success_codes = set(ordered_success_codes)
    failure_codes = set([FAIL_QC, FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG])
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
    if not variant_id in flags:
        print("Error, variant", variant_id, "missing from flags")
    assert(variant_id in flags)
    variant_flagged = False
    if flags[variant_id]["Filters"] != "PASS":
        variant_flagged = True
    return(variant_flagged)


def analyze_variant(variant, coverage_v2, coverage_v3, coverage_v4,
                    flags_v2, flags_v3,flags_v4, debug=False):
    """
    Analyze a single variant, adding the output columns
    """
    # Initialize the output columns.  First, check the read depth.
    # to see if the variant could have been observed (exception: non-exomic
    # regions for an exomic dataset).  Check to see if the variant was a SNV.
    # If the variant is observable and a SNV, initialize to PM2_Supporting
    # If it's observable and not a SNV, initialze to NO_CODE_NON_SNV
    # If it's not observable, initialize to FAIL_QC
    # If not observable and not a SNV, initialize to NO_CODE_NON_SNV
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
    (observable_in_v4,
     read_depth_v4) = estimate_coverage(int(variant["Hg38_Start"]),
                                        int(variant["Hg38_End"]),
                                        int(variant["Chr"]),coverage_v4,
                                        debug=debug)
    is_snv = (variant["Hg38_Start"] == variant["Hg38_End"]
              and len(variant["Ref"]) == 1 and len(variant["Alt"]) == 1)
    if is_snv:
        if observable_in_v2:
            variant_v2_code_id = PM2_SUPPORTING
        else:
            variant_v2_code_id = FAIL_QC
        if observable_in_v3:
            variant_v3_code_id = PM2_SUPPORTING
        else:
            variant_v3_code_id = NO_CODE
        if observable_in_v4:
            variant_v4_code_id = PM2_SUPPORTING
        else:
            variant_v4_code_id = NO_CODE
    else:
        variant_v2_code_id = NO_CODE_NON_SNV
        variant_v3_code_id = NO_CODE_NON_SNV
        variant_v4_code_id = NO_CODE_NON_SNV
    
    if debug:
        print("variant", variant["pyhgvs_cDNA"], " is snv:", is_snv, "preliminary codes", variant_v2_code_id,
              variant_v3_code_id, "observable", observable_in_v2, observable_in_v3)
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
        variant_v2_code_id = analyze_one_dataset(variant["faf95_popmax_exome_GnomAD"],
                                                         variant["faf95_popmax_population_exome_GnomAD"],
                                                         variant["Allele_count_exome_GnomAD"],
                                                         is_snv, read_depth_v2,
                                                         variant_is_flagged(variant["Variant_id_GnomAD"],
                                                                            flags_v2), debug)
        if debug:
            print("gnomAD V2 variant", variant["Variant_id_GnomAD"],
                  "popmax", variant["faf95_popmax_exome_GnomAD"],
                  "allele count", variant["Allele_count_exome_GnomAD"],
                  "read depth", read_depth_v2, "snv", is_snv,
                  "V2 code", variant_v2_code_id)
    variant_in_v3 = False
    if (field_defined(variant["Variant_id_GnomADv3"]) and observable_in_v3):
        variant_in_v3 = True
        variant_v3_code_id = analyze_one_dataset(variant["faf95_popmax_genome_GnomADv3"],
                                                         variant["faf95_popmax_population_genome_GnomADv3"],
                                                         variant["Allele_count_genome_GnomADv3"],
                                                         is_snv, read_depth_v3,
                                                         variant_is_flagged(variant["Variant_id_GnomADv3"],
                                                                            flags_v3), debug)
        if debug:
            print("gnomAD V3 variant", variant["Variant_id_GnomADv3"],
                  "popmax", variant["faf95_popmax_genome_GnomADv3"],
                  "allele count", variant["Allele_count_genome_GnomADv3"],
                  "read depth", read_depth_v3, "snv", is_snv,
                  "V3 code", variant_v3_code_id)
    variant_in_v4 = False
    if (field_defined(variant["Variant_id_GnomADv4"])):
        variant_in_v4 = True
        variant_v4_code_id = analyze_one_dataset(variant["faf95_popmax_joint_GnomADv4"],
                                                         variant["faf95_popmax_population_joint_GnomADv4"],
                                                         variant["Allele_count_joint_GnomADv4"],
                                                         is_snv, read_depth_v4,
                                                         variant_is_flagged(variant["Variant_id_GnomADv4"],
                                                                            flags_v4), debug)
        if debug:
            print("gnomAD V4 variant", variant["Variant_id_GnomADv4"],
                  "popmax", variant["faf95_popmax_joint_GnomADv4"],
                  "allele count", variant["Allele_count_joint_GnomADv4"],
                  "read depth", read_depth_v4, "snv", is_snv,
                  "V4 code", variant_v4_code_id, "V4 code", variant_v4_code_id)
    (variant[POPFREQ_CODE_ID],
     variant[POPFREQ_CODE_DESCR]) = analyze_across_datasets(variant_v2_code_id,variant["faf95_popmax_exome_GnomAD"],
                                                            variant["faf95_popmax_population_exome_GnomAD"],
                                                            variant_in_v2, variant_v3_code_id, 
                                                            variant["faf95_popmax_genome_GnomADv3"],
                                                            variant["faf95_popmax_population_genome_GnomADv3"],
                                                            variant_in_v3, is_snv, debug)
    variant["V4_Popfreq"] = variant_v4_code_id
    if debug:
        print("variant", variant["pyhgvs_cDNA"], "consensus code:", variant[POPFREQ_CODE_ID], "msg",
              variant[POPFREQ_CODE_DESCR], "v2 code", variant_v2_code_id,
              "v3 code", variant_v3_code_id)
    return()



def main():
    args = parse_args()
    cov2 = read_coverage(args.data_dir + '/df_cov_v2.csv')
    cov3 = read_coverage(args.data_dir + '/df_cov_v3.csv')
    cov4 = read_coverage(args.data_dir + '/df_cov_v4.csv')
    flags_v2 = read_flags(csv.DictReader(open(args.data_dir + "/brca.gnomAD.2.1.1.hg19.flags.tsv"),
                                         delimiter = "\t"))
    flags_v3 = read_flags(csv.DictReader(open(args.data_dir + "/brca.gnomAD.3.1.1.hg38.flags.tsv"),
                                         delimiter = "\t"))
    flags_v4 = read_flags(csv.DictReader(open(args.data_dir + "/brca.gnomAD.4.1.hg38.flags.tsv"),
                                         delimiter = "\t"))
    input_file = csv.DictReader(open(args.input), delimiter = "\t")
    output_file = initialize_output_file(input_file, args.output)
    for variant in input_file:
        if args.debug:
            print("About to analyze", variant["pyhgvs_cDNA"])
        analyze_variant(variant, cov2, cov3, cov4, flags_v2, flags_v3, flags_v4,
                        debug=args.debug)
        output_file.writerow(variant)


if __name__ == "__main__":
    main()

