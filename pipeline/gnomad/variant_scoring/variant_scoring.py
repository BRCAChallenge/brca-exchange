#!/usr/bin/env python
# coding: utf-8

import itertools
import argparse
import csv
from collections import OrderedDict

import numpy as np
import pandas as pd

#import gnomad.variant_scoring.constants as cnts
#from common import config, hgvs_utils, variant_utils, utils
#from data_merging.brca_pseudonym_generator import _normalize_genomic_fnc

GNOMAD_V2_CODE_ID = "Provisional_code_GnomAD"
GNOMAD_V2_CODE_DESCR = "Provisional_code_description_GnomAD"
GNOMAD_V3_CODE_ID = "Provisional_code_GnomADv3"
GNOMAD_V3_CODE_DESCR = "Provisional_code_description_GnomADv3"
POPFREQ_CODE_ID = "Provisional_evidence_code_popfreq"

FAIL_INSUFFICIENT_READ_DEPTH = "INSUFFICIENT_READ_DEPTH"
FAIL_VCF_FILTER_FLAG = "VCF_FILTER_FLAG"
FAIL_BOTH = "FAIL_BOTH"
FAIL_CONTRADICTORY = "FAIL_CONTRADICTORY_RESULTS"
BA1 = "BA1"
BS1 = "BS1"
BS1_SUPPORTING = "BS1_SUPPORTING"
NO_CODE = "NO_CODE"
PM2_SUPPORTING_ABSENT = "PM2_SUPPORTING (ABSENT)"
PM2_SUPPORTING_NOT_OBSERVED = "PM2_SUPPORTING (NOT OBSERVED)"
NEEDS_REVIEW = "NEEDS_REVIEW"

READ_DEPTH_THRESHOLD_FREQUENT_VARIANT = 20
READ_DEPTH_THRESHOLD_RARE_VARIANT = 25


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="built_final.tsv",
                        help="Input file with variant data")
    parser.add_argument("-o", "--output", default="built_with_popfreq.tsv",
                        help="Output file with new columns")
    parser.add_argument("-d", "--data_dir", default="./processed_brca",
                        help="Directory with the processed files")
    parser.add_argument("-r", "--resource_dir", default="BRCAExchange/resources",
                        help="Pipeline resources")
    args = parser.parse_args()
    return(args)


def var_obj_to_name(v):
    return '-'.join(str(x) for x in [v.chr, v.pos, v.ref, v.alt])

def var_name_to_obj(v):
    a = v.split('-')
    return variant_utils.VCFVariant(a[0], int(a[1]), a[2], a[3])


csep = '-'
def add_name_col(dfx: pd.DataFrame):
    dfx['var_name'] = (dfx['contigName'].apply(str) + csep +
                       (dfx['start']).apply(str) +
                       csep + dfx['referenceAllele'] + csep + dfx['alternateAlleles'])
    return dfx

def estimate_coverage(variant, start, end, chrom, df_cov):
    start = int(variant["pyhgvs_Hg38_Start"]) - 1
    end = int(variant["pyhgvs_Hg38_End"])
    positions = list(range(start, end))
    coverage_this_chrom = df_cov.loc[df_cov["chrom"] == int(variant["Chr"])]
    positions_this_variant = coverage_this_chrom[coverage_this_chrom["pos"].isin(positions)]
    meanval = positions_this_variant["mean"].mean()
    medianval = positions_this_variant["median"].median()
    return(meanval, medianval)


def coverage_per_variant(df_var: pd.DataFrame, df_cov: pd.DataFrame):
    """
    Computing median and mean coverage statistics per variant

    :param df_var: variant dataframe (from VCF)
    :param df_cov: dataframe with read coverage for every position
    :return:
    """
    df_var['positions'] = df_var.apply(lambda r: list(range(r['start'], r['end'])) , axis=1)
    df_var_with_read_stats = df_var.explode('positions').merge(df_cov, left_on=['contigName', 'positions'],
                                   right_on=['chrom', 'pos'] , how='left')
    cov_stats_per_var = (df_var_with_read_stats.rename(columns={'mean': 'pos_mean', 'median': 'pos_median'}).
            groupby('var_name').
            agg({'pos_mean': np.mean, 'pos_median': np.median}) # taking means of means and medians of medians for structural variants
            )
    df_var_with_cov = df_var.merge(cov_stats_per_var[['pos_mean', 'pos_median']], left_on='var_name', right_index=True)
    return df_var_with_cov


def calculate_faf95_cols(df_all):
    df_ret = df_all.copy()

    df_faf95 = df_ret[ cnts.faf95_col_names ]

    df_ret['faf95_popmax'] = df_faf95.idxmax(axis=1).apply(lambda s: s.split('_')[1] if s == s else s)
    df_ret['faf95_max'] = df_faf95.max(axis=1)

    return df_ret


def flatten_singleton_array_cols(df, col_names):
    df = df.copy()

    # precondition
    for c in col_names:
        assert df[c].apply(lambda x: len(x) == 1 if x is not None and not isinstance(x, float) else None).all(), f"column {c} does not contain singletons only"

    for c in col_names:
        df[c] = df[c].apply(lambda a: a[0] if a else None)

    return df


def aggregate_data(df_var, df_cov):
    # columns which have an 'array' type in the data frame, but actually contain just one entry
    singleton_cols = list(set(['alternateAlleles', 'popmax', 'filters', 'names', 'AN_popmax'] +
                          cnts.faf95_col_names +
                          [f"{pre}_{pop.lower()}" for pre in ['AC', 'AF'] for pop in cnts.pop_names + ['popmax']]
                          ))

    df_var = (flatten_singleton_array_cols(df_var,
                                           singleton_cols).
              rename(columns={'alternateAlleles ': 'alternateAllele'}))

    df_var = add_name_col(df_var)

    df_all = coverage_per_variant(df_var, df_cov)

    df_all = calculate_faf95_cols(df_all)

    df_all['max_faf_or_af'] = df_all['faf95_max'].fillna(df_all['AF_popmax'])

    df_all['max_faf_or_af_pop'] = df_all['faf95_popmax'].fillna(df_all['popmax'])

    df_all['is_snv'] = (df_all['end'] - df_all['start']) == 1

    return df_all


def add_filter_flag_col(df, cols):
    df['vcf_filter_flag'] = df['filters'].apply(
        lambda filt_el: len(set(filt_el.split(';')).intersection(cols)) > 0)


def add_read_depth_col(df, threshold):
    df['sufficient_read_depth'] = (df['pos_mean'] >= threshold) & (
            df['pos_median'] >= threshold)


def process_v2(df_var2, df_cov2, read_depth_thresh, resource_dir):
    agg_v2 = aggregate_data(df_var2, df_cov2)

    agg_v2['src'] = 'v2'
    add_filter_flag_col(agg_v2, set(['InbreedingCoeff', 'RF', 'AC0']))
    add_read_depth_col(agg_v2, read_depth_thresh)

    variants = [var_name_to_obj(v) for v in agg_v2['var_name']]

    # assuming same layout of resource_dir as for the main pipeline
    chain_file = resource_dir / 'hg19ToHg38.over.chain.gz'
    ref_file = resource_dir / 'hg38.fa'
    vars_gnomad2_hg38, failed_entries = variant_utils.convert_to_hg38(variants, chain_file, ref_file, resource_dir)

    assert len(agg_v2) == len(vars_gnomad2_hg38) + len(failed_entries)

    failed_keys = [var_obj_to_name(v) for v in failed_entries]
    failed_indices = set(agg_v2[agg_v2['var_name'].isin(failed_keys)].index)
    all_indices = set(agg_v2.index)

    valid_indices = all_indices - failed_indices

    entries_sorted = sorted(itertools.chain(zip(sorted(valid_indices), vars_gnomad2_hg38),
                                            zip(failed_indices, [None] * len(failed_indices))), key=lambda p: p[0])

    agg_v2_lift_over = agg_v2.rename(columns={'var_name': 'var_name_hg19'}).assign(
        var_name=lambda _: pd.Series(var_obj_to_name(p[1]) if p[1] else p[1] for p in entries_sorted))

    return agg_v2_lift_over


def process_v3(df_var3, df_cov3, read_depth_thresh):
    agg_v3 = aggregate_data(df_var3, df_cov3)

    agg_v3['src'] = 'v3'
    add_filter_flag_col(agg_v3, set(['AS_VQSR', 'AC0']))
    add_read_depth_col(agg_v3, read_depth_thresh)
    return agg_v3


def determine_evidence_code_per_variant(r):
    if not r['sufficient_read_depth']:
        return 'fail_insufficient_read_depth'

    if r['vcf_filter_flag']:
        return 'fail_vcf_filter_flag'

    if r['max_faf_or_af'] > 0.001:
        return 'BA1'

    if r['max_faf_or_af'] > 0.0001:
        return 'BS1'

    if not np.isnan(r['max_faf_or_af']):
        # low max_faf_or_af
        return 'code_missing'

    if r['is_snv']:
        return 'pm2_supporting'

    return 'need_review'


def add_final_code_column(df):
    success_codes = set(['BA1', 'BS1', 'pm2_supporting'])
    dq_failure_codes = set(['fail_insufficient_read_depth', 'fail_vcf_filter_flag', 'need_review'])
    CODE_MISSING = 'code_missing'
    
    def is_both_sources(dfg):
        return len(dfg) > 1

    def set_final_code(dfg):
        if len(dfg) == 1:
            return dfg['evidence_code'].iloc[0]
        elif len(dfg) == 2:
            # we have data from both v2 and v3
            c1 = dfg['evidence_code'].iloc[0]
            c2 = dfg['evidence_code'].iloc[1]

            if (c1 in dq_failure_codes and c2 == CODE_MISSING) or (c2 in dq_failure_codes and c1 == CODE_MISSING):
                return CODE_MISSING
            
            if c1 in dq_failure_codes and c2 in dq_failure_codes:
                return "fail_both"
            
            if (c1 == 'pm2_supporting' and c2 == CODE_MISSING) or (c2 == 'pm2_supporting' and c1 == CODE_MISSING):
                return CODE_MISSING

            if c1 == CODE_MISSING and c2 == CODE_MISSING:
                return CODE_MISSING
            
            if c2 in success_codes and c1 not in success_codes:
                return c2
            if c1 in success_codes and c2 not in success_codes:
                return c1

            assert c1 in success_codes and c2 in success_codes
            
            if c1 == c2:
                return c1
            else:
                return "contradictory_results"

        raise ValueError("some duplicate variants per source?")

    v2_and_v3 = df.groupby('var_name').apply(is_both_sources)
    per_variant_code = df.groupby('var_name').apply(set_final_code)

    return df.merge(pd.DataFrame({'in_v2_and_v3' : v2_and_v3, 'final_code': per_variant_code}).reset_index(), how='left')


def extract_variant_scoring_data(df_cov2, df_cov3, df_var2, df_var3, read_depth_thresh, resource_dir):
    agg_v2 = process_v2(df_var2, df_cov2, read_depth_thresh, resource_dir)
    agg_v3 = process_v3(df_var3, df_cov3, read_depth_thresh)

    df_overall = pd.concat([agg_v2, agg_v3]).reset_index(drop=True)

    # running the actual scoring algorithm and the combined data
    df_overall['evidence_code'] = df_overall.apply(determine_evidence_code_per_variant, axis=1)
    df_overall = add_final_code_column(df_overall)

    return df_overall.sort_values('var_name')


def add_normalization_cols(df, strand_dict, processes=2):
    TMP_HGVS_HG38 = 'tmp_hgvs_hg38'
    TMP_VAR_OBJ_FIELD = 'tmp_var_obj'
    TMP_GENE_SYMBOL = 'Gene_Symbol'

    hgvs_proc = hgvs_utils.HgvsWrapper()

    df[TMP_GENE_SYMBOL] = df['contigName']
    df[TMP_VAR_OBJ_FIELD] = df['var_name'].apply(lambda v: var_name_to_obj(v))
    df[TMP_HGVS_HG38] = df[TMP_VAR_OBJ_FIELD].apply(
        lambda v: v.to_hgvs_obj(hgvs_proc.contig_maps[hgvs_utils.HgvsWrapper.GRCh38_Assem]))

    df = utils.parallelize_dataframe(df, _normalize_genomic_fnc(TMP_HGVS_HG38,
                                                                'var_name_right', True,
                                                                strand_dict), processes)
    df = utils.parallelize_dataframe(df, _normalize_genomic_fnc(TMP_HGVS_HG38,
                                                                'var_name_left', False,
                                                                strand_dict), processes)

    def _convert_to_name(converted):
        return converted.apply(lambda h: var_obj_to_name(variant_utils.VCFVariant.from_hgvs_obj(h)))

    df['var_name_right'] = _convert_to_name(df['var_name_right'])
    df['var_name_left'] = _convert_to_name(df['var_name_left'])

    return df.drop(columns=[TMP_HGVS_HG38, TMP_VAR_OBJ_FIELD, TMP_GENE_SYMBOL])

def initialize_output_file(input_file, output_filename):
    """
    Create an empty output file with the new columns
    """
    new_columns = [GNOMAD_V2_CODE_ID, GNOMAD_V2_CODE_DESCR,
                   GNOMAD_V3_CODE_ID, GNOMAD_V3_CODE_DESCR,
                   POPFREQ_CODE_ID]
    input_header_row = input_file.fieldnames
    output_header_row = input_header_row + new_columns
    output_file = csv.DictWriter(open(output_filename,"w"), fieldnames=output_header_row,
                                 delimiter = '\t')
    output_file.writeheader()
    return(output_file)

def field_defined(field):
    """
    Return a binary value indicating whether or not this variant has the popmax FAF defined
    """
    return(field != "-")

def analyze_dataset(faf95_popmax_str, faf95_population, allele_count, is_snv, mean_read_depth, median_read_depth, debug=True):
    #
    # Get the coverage data
    #if not r['sufficient_read_depth']:
    #    return 'fail_insufficient_read_depth'
    #
    read_depth = min(mean_read_depth, median_read_depth)
    #
    # Get the VCF filter flags somehow
    #if r['vcf_filter_flag']:
    #    return 'fail_vcf_filter_flag'
    #
    vcf_filter_flag = False  ### HERE
    #
    rare_variant = False
    if field_defined(faf95_popmax_str):
        faf = float(faf95_popmax_str)
        if np.isnan(faf):
            rare_variant = True
        elif faf <= 0.00002:
            rare_variant = True
    else:
        faf = None
        rare_variant = True
    if rare_variant and read_depth < READ_DEPTH_THRESHOLD_RARE_VARIANT:
        return(FAIL_INSUFFICIENT_READ_DEPTH, "Insufficient read depth")
    if (not rare_variant) and read_depth < READ_DEPTH_THRESHOLD_FREQUENT_VARIANT:
        return(FAIL_INSUFFICIENT_READ_DEPTH, "Insufficient read depth")
    #
    # Address the cases where a variant cannot be analyzed because the gnomAD data is flagged
    if vcf_filter_flag:
        return(FAIL_VCF_FILTER_FLAG, "Data flagged in the gnomAD VCF")
    #
    # Address the cases where FAF is defined, and the variant is a candidate for a
    # evidence code for high population frequency (BA1, BS1, BS1_SUPPORTING)
    if not rare_variant:
        if faf > 0.001:
            return(BA1, "FAF > 0.001 in " + faf95_population)
        elif faf >  0.0001:
            return(BS1, "FAF > 0.0001 in " + faf95_population)
        elif faf > 0.00002:
            return(BS1_SUPPORTING, "FAF > 0.00002 in " + faf95_population)
    if rare_variant:
        if not field_defined(allele_count):
            return(NO_CODE, "Variant is rare in the dataset")
        elif int(allele_count) > 0:
            return(NO_CODE, "Variant is rare in the dataset")
        else:
            if is_snv:
                return(PM2_SUPPORTING_NOT_OBSERVED, "Variant absent from the dataset")
            else:
                return(NO_CODE, "Rare indel or structural variant")
    return(NEEDS_REVIEW, "This variant_needs_review")


def analyze_across_datasets(code_v2, code_v3):
    """
    Given the per-dataset evidence codes, generate an overall evidence code
    """
    benign_codes = [BA1, BS1, BS1_SUPPORTING]
    pathogenic_codes = [PM2_SUPPORTING_NOT_OBSERVED]
    intermediate_codes = [ NO_CODE]
    ordered_success_codes = benign_codes + intermediate_codes + pathogenic_codes
    success_codes = set(ordered_success_codes)
    failure_codes = set([FAIL_INSUFFICIENT_READ_DEPTH, FAIL_VCF_FILTER_FLAG, NEEDS_REVIEW])
    #
    # If the variant is absent in one of the two datasets, report its code from the one where
    # it is present.  If the variant is absent in both datasets, the return code will be
    # PM2_SUPPORTING_ABSENT
    if code_v2 == PM2_SUPPORTING_ABSENT and code_v3 != PM2_SUPPORTING_ABSENT:
        return(code_v3)
    elif code_v3 == PM2_SUPPORTING_ABSENT:
        #
        # Here, either code_v2 is PM2_SUPPORTING_ABSENT or it's not.  Either way, the correct
        # action is to return code_v2.  If it's PM2_SUPPORTING_ABSENT, then it's PM2_SUPPORTING_ABSENT
        # in both.  If it's not PM2_SUPPORTING_ABSENT, then the correct output is whatever code_v2 is.
        return(code_v2)
    #
    # If the variant had an error in both datasets, report that it had an error.  If the variant
    # had an error in only one of the two datwsets, report the one in which it could be called.
    if code_v2 in failure_codes and code_v3 in failure_codes:
        return(FAIL_BOTH)
    elif code_v2 in failure_codes:
        return(code_v3)
    elif code_v3 in failure_codes:
        return(code_v2)
    #
    # At this point, it should be true that the variant has been observed in both datasets
    # and neither one had a failure code.
    assert(code_v2 in ordered_success_codes)
    assert(code_v3 in ordered_success_codes)
    #
    # If the variant had a pathogenic code in one dataset and a benign code in the other,
    # flag it for review
    if code_v2 in benign_codes and code_v3 in pathogenic_codes:
        return(FAIL_CONTRADICTORY)
    elif code_v2 in pathogenic_codes and code_v3 in benign_codes:
        return(FAIL_CONTRADICTORY)
    #
    # If the codes are not contradictory, then return the stronger result
    if ordered_success_codes.index(code_v2) < ordered_success_codes.index(code_v3):
        return(code_v2)
    else:
        return(code_v3)


def analyze_variant(variant, coverage_v2, coverage_v3, debug=True):
    """
    Analyze a single variant, adding the output columns
    """
    # Initialize the output columns with the assumption that the variant isn't observed
    variant[GNOMAD_V2_CODE_ID] = PM2_SUPPORTING_ABSENT
    variant[GNOMAD_V2_CODE_DESCR] = "Variant not observed in the gnomAD V2 Exomes"
    variant[GNOMAD_V3_CODE_ID] = PM2_SUPPORTING_ABSENT
    variant[GNOMAD_V3_CODE_DESCR] = "Variant not observed in the gnomAD V3 Genomes"
    variant[POPFREQ_CODE_ID] = PM2_SUPPORTING_ABSENT
    is_snv = (variant["Hg38_Start"] == variant["Hg38_End"])
    if debug:
        print("variant is snv:", is_snv)
    if field_defined(variant["Variant_id_GnomAD"]):
        (mean_read_depth, median_read_depth) = estimate_coverage(int(variant["pyhgvs_Hg19_Start"]),
                                                                 int(variant["pyhgvs_Hg19_End"]),
                                                                 int(variant["Chr"]),coverage_v2)
        if debug:
            print("gnomAD V2 variant", variant["Variant_id_GnomAD"], "popmax", variant["faf95_popmax_exome_GnomAD"],
                  "allele count", variant["Allele_count_exome_GnomAD"], "mean read depth", mean_read_depth,
                  "median read depth", median_read_depth)
        (variant[GNOMAD_V2_CODE_ID],
         variant[GNOMAD_V2_CODE_DESCR]) = analyze_dataset(variant["faf95_popmax_exome_GnomAD"],
                                                          variant["faf95_popmax_population_exome_GnomAD"],
                                                          variant["Allele_count_exome_GnomAD"],
                                                          is_snv, mean_read_depth, median_read_depth)
        if debug:
            print("From gnomAD V2: code ID", variant[GNOMAD_V2_CODE_ID], "descr", variant[GNOMAD_V2_CODE_DESCR])
    if field_defined(variant["Variant_id_GnomADv3"]):
        (mean_read_depth, median_read_depth) = estimate_coverage(int(variant["pyhgvs_Hg38_Start"]),
                                                                 int(variant["pyhgvs_Hg38_End"]),
                                                                 int(variant["Chr"]),coverage_v3)
        if debug:
            print("gnomAD V3 variant", variant["Variant_id_GnomADv3"], "popmax", variant["faf95_popmax_genome_GnomADv3"],
                  "allele count", variant["Allele_count_genome_GnomADv3"], "mean read depth", mean_read_depth,
                  "median read depth", median_read_depth)
        (variant[GNOMAD_V3_CODE_ID],
         variant[GNOMAD_V3_CODE_DESCR]) = analyze_dataset(variant["faf95_popmax_genome_GnomADv3"],
                                                          variant["faf95_popmax_population_genome_GnomADv3"],
                                                          variant["Allele_count_genome_GnomADv3"],
                                                          is_snv, mean_read_depth, median_read_depth)
        if debug:
            print("From gnomAD V3: code ID", variant[GNOMAD_V3_CODE_ID], "descr", variant[GNOMAD_V3_CODE_DESCR])
    variant[POPFREQ_CODE_ID] = analyze_across_datasets(variant[GNOMAD_V2_CODE_ID],
                                                        variant[GNOMAD_V3_CODE_ID])
    if debug:
        print("consensus code:", variant[POPFREQ_CODE_ID])
    return()


def main():
    args = parse_args()
    print(args)
    #cfg_df = config.load_config(gene_config_path)
    df_cov2 = pd.read_parquet(args.data_dir + '/df_cov_v2.parquet')
    df_cov3 = pd.read_parquet(args.data_dir + '/df_cov_v3.parquet')
    input_file = csv.DictReader(open(args.input), delimiter = "\t")
    output_file = initialize_output_file(input_file, args.output)
    for variant in input_file:
        analyze_variant(variant, df_cov2, df_cov3)
        output_file.writerow(variant)
                
    #df_var2 = pd.read_parquet(data_dir / 'df_var_v2.parquet')
    #df_var3 = pd.read_parquet(data_dir / 'df_var_v3.parquet')

    #read_depth_thresh = 30

    #df = extract_variant_scoring_data(df_cov2, df_cov3, df_var2, df_var3, read_depth_thresh, Path(resource_dir))

    # add var_name columns with different normalization to join data with other sources (e.g. brca exchange output data)
    #strand_dict = { int(r['chr']) : r[config.STRAND_COL] for _, r in cfg_df.iterrows() }
    #df = add_normalization_cols(df, strand_dict)

    #df.to_parquet(output_path)


if __name__ == "__main__":
    main()

