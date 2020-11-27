import importlib
import logging
from collections import namedtuple
from typing import List, Callable, Sequence, TypeVar

import numpy as np
import pandas as pd
import pathos
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
                           chr_intervals).items():
        interval_tuples = [
            (r.start, r.end + 1, interval_tuple_builder(c, r.start, r.end)) for
            r in regs]

        d[c] = IntervalTree.from_tuples(
            interval_tuples)

    return d


T = TypeVar('T')


def split_list_in_chunks(seq: Sequence[T], chunk_size: int) -> List[List[T]]:
    return [seq[pos:pos + chunk_size] for pos in range(0, len(seq), chunk_size)]


def parallelize_dataframe(df: pd.DataFrame, func: Callable[[pd.DataFrame], pd.DataFrame], n_processes: int = 4):
    if df.empty:
        return df

    df_split = np.array_split(df, n_processes)

    with pathos.multiprocessing.ProcessingPool(ncpus=n_processes) as pool:
        df = pd.concat(pool.map(func, df_split))

    return df


def setup_logfile(log_path, log_level=logging.INFO):
    # https://stackoverflow.com/questions/20240464/python-logging-file-is-not-working-when-using-logging-basicconfig
    importlib.reload(logging)

    logging.basicConfig(filename=log_path, filemode="w", level=log_level,
                        format=' %(asctime)s %(filename)-15s %(message)s')
