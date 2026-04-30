#!/usr/bin/env python
# coding: utf-8

import argparse
import bisect
import csv
import math
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

 
NO_CODE_NON_SNV_MSG = (f"This [duplication/insertion/deletion/delins/large genomic rearrangement] variant "
                       f"was not observed in gnomAD v4.1, but PM2_Supporting was not applied since recall "
                       f"is considered suboptimal for this type of variant (PM2_Supporting not met). ")
PM2_SUPPORTING_MSG = "This variant is absent from gnomAD v4.1 (PM2_Supporting met)."
FAIL_INSUFFICIENT_READ_DEPTH_OR_FILTER_FLAG_MSG = (f"This variant is present in gnomAD but is not meeting the "
                                                    f"specified read depths threshold \u2265{READ_DEPTH_THRESHOLD_FREQUENT_VARIANT} "
                                                    f"(PM2_Supporting, BS1, and BA1 are not met).")
FAIL_LCR_MSG = "PM2_Supporting was not applied since this variant overlaps a low-complexity region (LCR)."

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="built_final.tsv",
                        help="Input file with variant data")
    parser.add_argument("-o", "--output", default="built_with_popfreq.tsv",
                        help="Output file with new columns")
    parser.add_argument("-d", "--data_dir", default="./processed_brca",
                        help="Directory with the processed files")
    parser.add_argument("--popfreq-code-id", default=POPFREQ_CODE_ID,
                        help="Column header for the v4 popfreq evidence code")
    parser.add_argument("--popfreq-code-descr", default=POPFREQ_CODE_DESCR,
                        help="Column header for the v4 popfreq evidence description")
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




def initialize_output_file(input_file, output_filename,
                           popfreq_code_id, popfreq_code_descr):
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

    

def analyze_one_dataset(faf95_popmax_str, allele_count, snv_or_small_indel,
                        read_depth, vcf_filter_flag,
                        allele_count_threshold,
                        chrom, genome_start, genome_end, lcr,
                        faf95_popmax, faf95_popmax_population,
                        allele_count_pop, allele_number_pop, debug=True):
    print("debug", debug)
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


def get_population_allele_counts(population, variant):
    """
    Given a population code (e.g. "AMR") from the gnomAD v4 joint data,
    return the population-specific allele count and allele number.
    Returns (None, None) if population is None or "NA", or if the
    expected keys are absent from the variant dictionary.
    """
    if population is None or population == "NA":
        return None, None
    s1 = f"Allele_count_joint_{population}_GnomADv4"
    s2 = f"Allele_number_joint_{population}_GnomADv4"
    if s1 in variant and s2 in variant:
        return variant[s1], variant[s2]
    return None, None


def variant_is_flagged(variant_id, flags):
    if not field_defined(variant_id):
        return False
    if not variant_id in flags:
        print("Error, variant", variant_id, "missing from flags")
    assert(variant_id in flags)
    variant_flagged = False
    if flags[variant_id]["Filters"] != "PASS":
        variant_flagged = True
    return(variant_flagged)


def analyze_variant(variant, coverage_v4, flags_v4,
                    popfreq_code_id=POPFREQ_CODE_ID,
                    popfreq_code_descr=POPFREQ_CODE_DESCR,
                    lcr=None, debug=False):
    """
    Analyze a single variant, adding the output columns
    """
    (observable_in_v4,
     read_depth_v4) = estimate_coverage(int(variant["Hg38_Start"]),
                                        int(variant["Hg38_End"]),
                                        int(variant["Chr"]), coverage_v4,
                                        debug=debug)
    snv_or_small_indel = (int(variant["Hg38_End"]) - int(variant["Hg38_Start"]) <= SMALL_INDEL_SIZE_THRESHOLD)

    if debug:
        print("variant", variant["pyhgvs_cDNA"], " snv or small indel:", snv_or_small_indel)

    variant[popfreq_code_id] = ""
    variant[popfreq_code_descr] = ""

    population = variant["faf95_popmax_population_joint_GnomADv4"]
    allele_count_pop, allele_number_pop = get_population_allele_counts(population, variant)
    (code_id, msg) = analyze_one_dataset(variant["faf95_popmax_joint_GnomADv4"],
                                         variant["Allele_count_joint_GnomADv4"],
                                         snv_or_small_indel, read_depth_v4,
                                         variant_is_flagged(variant["Variant_id_GnomADv4"], flags_v4),
                                         ALLELE_COUNT_RARE_VARIANT_THRESHOLD,
                                         int(variant["Chr"]), int(variant["Hg38_Start"]),
                                         int(variant["Hg38_End"]), lcr,
                                         variant["faf95_popmax_joint_GnomADv4"],
                                         variant["faf95_popmax_population_joint_GnomADv4"],
                                         allele_count_pop, allele_number_pop,
                                         debug=debug)
    if debug:
        print("gnomAD V4 variant", variant["Variant_id_GnomADv4"],
              "popmax", variant["faf95_popmax_joint_GnomADv4"],
              "allele count", variant["Allele_count_joint_GnomADv4"],
              "read depth", read_depth_v4, "snv", snv_or_small_indel,
              "flagged", variant_is_flagged(variant["Variant_id_GnomADv4"], flags_v4),
              "V4 code", code_id)
    variant[popfreq_code_id] = code_id
    variant[popfreq_code_descr] = msg
    if debug:
        print("variant", variant["pyhgvs_cDNA"], "v4 code:", variant[popfreq_code_id])
    return()



def main():
    args = parse_args()
    cov4 = read_coverage(args.data_dir + '/df_cov_v4.csv')
    flags_v4 = read_flags(csv.DictReader(open(args.data_dir + "/brca.gnomAD.4.1.hg38.flags.tsv"),
                                         delimiter = "\t"))
    lcr = read_lcr(args.lcr) if args.lcr else None
    input_file = csv.DictReader(open(args.input), delimiter = "\t")
    output_file = initialize_output_file(input_file, args.output,
                                         args.popfreq_code_id,
                                         args.popfreq_code_descr)
    for variant in input_file:
        if args.debug:
            print("About to analyze", variant["pyhgvs_cDNA"])
        analyze_variant(variant, cov4, flags_v4,
                        popfreq_code_id=args.popfreq_code_id,
                        popfreq_code_descr=args.popfreq_code_descr,
                        lcr=lcr, debug=args.debug)
        output_file.writerow(variant)


if __name__ == "__main__":
    main()

