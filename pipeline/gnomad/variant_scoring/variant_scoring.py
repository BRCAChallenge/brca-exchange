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

import gnomad.variant_scoring.constants as cnts

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
    df_faf95 = df_all[cnts.faf95_col_names]
    df_all['faf95_popmax'] = df_faf95.idxmax(axis=1).apply(lambda s: s.split('_')[1] if s == s else s)
    df_all['faf95_max'] = df_faf95.max(axis=1)


def aggregate_data(df_var, df_cov):
    df_all = coverage_per_variant(df_var, df_cov)

    # TODO: use return!
    calculate_faf95_cols(df_all)

    df_all['max_faf_or_af'] = df_all['faf95_max'].fillna(df_all['AF_popmax'])

    df_all['max_faf_or_af_pop'] = df_all['faf95_popmax'].fillna(df_all['popmax'])

    return df_all
