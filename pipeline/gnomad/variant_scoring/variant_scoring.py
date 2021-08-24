#!/usr/bin/env python
# coding: utf-8

import itertools
from pathlib import Path

import click
import numpy as np
import pandas as pd

import gnomad.variant_scoring.constants as cnts
from common import variant_utils


def var_obj_to_name(v):
    return '-'.join(str(x) for x in [v.chr, v.pos, v.ref, v.alt])

def var_name_to_obj(v):
    a = v.split('-')
    return variant_utils.VCFVariant(a[0], int(a[1]), a[2], a[3])


csep = '-'
def add_name_col(dfx: pd.DataFrame):
    dfx['var_name'] = (dfx['contigName'].apply(str) + csep +
                       (dfx['start'] + 1).apply(str) +  # used rs* variant identifiers at https://gnomad.broadinstitute.org/ to verify names correspond
                       csep + dfx['referenceAllele'] + csep + dfx['alternateAlleles'])
    return dfx


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

    # TODO: move out and generalize
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


def determine_evidence_code(r):
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
        return 'pm2_supporting'  # TODO: correct?

    return 'need_review'


def extract_variant_scoring_data(df_cov2, df_cov3, df_var2, df_var3, read_depth_thresh, resource_dir):
    agg_v2 = process_v2(df_var2, df_cov2, read_depth_thresh, resource_dir)
    agg_v3 = process_v3(df_var3, df_cov3, read_depth_thresh)

    df_overall = pd.concat([agg_v2, agg_v3]).sort_values('var_name').reset_index(drop=True)

    df_overall['evidence_code'] = df_overall.apply(determine_evidence_code, axis=1)

    return df_overall


@click.command()
@click.argument('data_dir', type=click.Path(readable=True))
@click.argument('output_path', type=click.Path(writable=True))
@click.option('--resource-dir', type=click.Path(readable=True), help="resource dir for lift over (same as for the main pipeline)")
def main(data_dir, output_path, resource_dir):
    data_dir = Path(data_dir)
    df_cov2 = pd.read_parquet(data_dir / 'df_cov_v2.parquet')
    df_cov3 = pd.read_parquet(data_dir / 'df_cov_v3.parquet')
    df_var2 = pd.read_parquet(data_dir / 'df_var_v2.parquet')
    df_var3 = pd.read_parquet(data_dir / 'df_var_v3.parquet')

    read_depth_thresh = 30 # TODO parametrize

    df = extract_variant_scoring_data(df_cov2, df_cov3, df_var2, df_var3, read_depth_thresh, Path(resource_dir))


    df.to_parquet(output_path)


if __name__ == "__main__":
    main()

