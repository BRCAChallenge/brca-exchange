import pandas as pd

from math import floor, log10
from collections import namedtuple

from intervaltree import IntervalTree
from toolz import groupby

ChrInterval = namedtuple("ChrInterval", "chr, start, end")


def build_interval_trees_by_chr(chr_intervals, interval_tuple_builder):
    '''
    Build a dictionary from chromosome to a interval tree.
    The interval tree may contain some additional payload which can be added a additional
    tuple element in the interval_tuple_builder

    :param chr_intervals: Iterable[ChrInterval]
    :param interval_tuple_builder: Function[ChrInterval, object]
    :return: dict[int, IntervalTree]
    '''
    d = {}

    for c, regs in groupby(lambda r: r.chr,
                           chr_intervals).iteritems():
        interval_tuples = [
            (r.start, r.end + 1, interval_tuple_builder(c, r.start, r.end)) for
            r in regs]

        d[c] = IntervalTree.from_tuples(
            interval_tuples)

    return d


def get_genome_regions_symbol_dict(gene_config_df):
    '''
    Build a dictionary from chromosome to a interval tree with the gene symbol
    as payload.

    Used to return the gene symbol (e.g. BRCA1) of a given position.

    :param gene_config_df:
    :return: dict[int, IntervalTree]
    '''

    def interval_tree_mapper(c, s, e):
        return gene_regions_dict[ChrInterval(c, s, e)]['symbol']

    gene_regions_dict = extract_gene_regions_dict(gene_config_df)

    return build_interval_trees_by_chr(gene_regions_dict.keys(),
                                    interval_tree_mapper)


def extract_gene_regions_dict(gene_config_df):
    '''
    Prepares dict used in get_genome_regions_symbol_dict

    :param gene_config_df: gene metadata dataframe
    :return: dict[ChrInterval, dict[str, object]]
    '''

    # a[2]+1 -> end_hg38 position is inclusive, but in seq_utils end position is treated exclusive
    return {ChrInterval(a[0], a[1], a[2] + 1): {'symbol': a[3]} for a
            in
            gene_config_df.loc[:,
            ['chr', 'start_hg38', 'end_hg38', 'symbol']].values}

def load_config(path):
    '''
    Loads gene metadata from file into a dataframe

    :param path: config file path
    :return: dataframe
    '''
    df = pd.read_csv(path, sep=',', header=0)
    return df.set_index('symbol', drop=False)


def round_sigfigs(num, sig_figs):
    if num != 0:
        return round(num, -int(floor(log10(abs(num))) - (sig_figs - 1)))
    else:
        return 0  # Can't take the log of 0
