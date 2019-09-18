from multiprocessing.pool import ThreadPool

import numpy as np
import pandas as pd


def parallelize_dataframe(df, func, n_cores=4):
    df_split = np.array_split(df, n_cores)

    pool = ThreadPool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df
