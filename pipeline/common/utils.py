from collections import namedtuple
from multiprocessing.pool import ThreadPool

import numpy as np
import pandas as pd
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


def parallelize_dataframe(df, func, n_cores=4):
    df_split = np.array_split(df, n_cores)

    pool = ThreadPool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df
