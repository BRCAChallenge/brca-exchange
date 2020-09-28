import pandas as pd

from common import utils


def test_parallelize_dataframe():
    def double_numbers(df):
        df['doubled'] = 2 * df['numbers']
        return df

    pd.testing.assert_frame_equal(utils.parallelize_dataframe(pd.DataFrame(), double_numbers), pd.DataFrame())

    nrows = 10
    df = pd.DataFrame({"numbers": range(0, nrows)})
    expected = 2 * df['numbers']

    pd.testing.assert_series_equal(utils.parallelize_dataframe(df, double_numbers)['doubled'], expected,
                                   check_names=False)
