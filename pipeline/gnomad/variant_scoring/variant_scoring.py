#!/usr/bin/env python
# coding: utf-8

import glob
import math
import numpy as np
import os
import pandas as pd
import re
import shutil
import sys
from pathlib import Path
import itertools

import gnomad.variant_scoring.constants as cnts
from common import variant_utils

resource_dir = Path('/Users/marc/nobackup/brca/wdir/resources')
chain_file = resource_dir/'hg19ToHg38.over.chain.gz'
ref_file = resource_dir/'hg38.fa'


def var_obj_to_name(v):
    return '-'.join(str(x) for x in [v.chr, v.pos, v.ref, v.alt])

def var_name_to_obj(v):
    a = v.split('-')
    return variant_utils.VCFVariant(a[0], a[1], a[2], a[3])

csep = '-'
def add_name_col(df_cov: pd.DataFrame):
    df_cov['var_name'] = (df_cov['contigName'].apply(str) + csep + df_cov['start'].apply(str) + csep
                       + df_cov['referenceAllele'] + csep + df_cov['alternateAlleles'])
    return df_cov


def coverage_per_variant(df_var: pd.DataFrame, df_cov: pd.DataFrame):
    df_var['intervals'] = df_var.apply(lambda r: list(range(r['start'], r['end'])) , axis=1)

    df_ex = df_var.explode('intervals').merge(df_cov, left_on=['contigName', 'intervals'],
                                   right_on=['chrom', 'pos'] , how='left')

    dfx =(df_ex.rename(columns = {'mean':'pos_mean', 'median' : 'pos_median' }).
            groupby('var_name')[['pos_mean', 'pos_median']].agg(['mean', 'median'])
          )

    dfx.columns.set_levels([f"{l1}_{l2}" for l1, l2 in zip(dfx.columns.levels[0], dfx.columns.levels[1])], level=0)


    dfxs = dfx[[('pos_mean', 'mean'), ('pos_median', 'median')]]
    dfxs.columns = dfxs.columns.droplevel(1)

    df_var_with_cov = df_var.merge(dfxs[['pos_mean', 'pos_median']], left_on='var_name', right_index=True)

    return df_var_with_cov


def calculate_faf95_cols(df_all):
    df_faf95 = df_all[ cnts.faf95_col_names ]
    df_all['faf95_popmax'] = df_faf95.idxmax(axis=1).apply(lambda s: s.split('_')[1] if s == s else s)
    df_all['faf95_max'] = df_faf95.max(axis=1)

def extract_singleton_array_cols(df, col_names):
    # precondition
    for c in col_names:
        assert df[c].apply(lambda x: len(x) == 1 if x else None).all(), f"column {c} does not contain singletons only"

    for c in col_names:
        df[c] = df[c].apply(lambda a: a[0] if a else None)

    return df


def aggregate_data(df_var, df_cov):
    df_var = (extract_singleton_array_cols(df_var, ['alternateAlleles', 'popmax'] + cnts.faf95_col_names).
              rename(columns={'alternateAlleles ': 'alternateAllele'}))

    df_var = add_name_col(df_var)


    df_all = coverage_per_variant(df_var, df_cov)

    # TODO: use return!
    calculate_faf95_cols(df_all)

    df_all['max_faf_or_af'] = df_all['faf95_max'].fillna(df_all['AF_popmax'])

    df_all['max_faf_or_af_pop'] = df_all['faf95_popmax'].fillna(df_all['popmax'])

    return df_all



def add_filter_flag_col(df, cols):
    df['vcf_filter_flag'] = df['filters'].apply(
        lambda filt_el: len(set(filt_el[0].split(';')).intersection(cols)) > 0)

def add_read_depth_col(df, threshold):
    df['sufficient_read_depth'] = (df['pos_mean'] >= threshold) & (
            df['pos_median'] >= threshold)

def process_v2(df_var2, df_cov2):
    read_depth_thresh = 30 # TODO: parametrize
    agg_v2 = aggregate_data(df_var2, df_cov2)

    # crossmap!

    agg_v2['src'] = 'v2'
    add_filter_flag_col(agg_v2, set(['InbreedingCoeff', 'RF', 'AC0']))
    add_read_depth_col(agg_v2, read_depth_thresh)

    variants = [var_name_to_obj(v) for v in agg_v2['var_name']]

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


def process_v3(df_var3, df_cov3):
    read_depth_thresh = 30  # TODO: parametrize
    agg_v3 = aggregate_data(df_var3, df_cov3)

    # TODO: what about InbreedingCoeff (2 cases)?
    agg_v3['src'] = 'v3'
    add_filter_flag_col(agg_v3, set(['AS_VQSR', 'AC0']))
    add_read_depth_col(agg_v3, read_depth_thresh)
    return agg_v3
