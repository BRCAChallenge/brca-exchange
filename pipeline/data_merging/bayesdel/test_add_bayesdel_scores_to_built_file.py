import os
import pandas as pd
import numpy as np
from click.testing import CliRunner

from data_merging.bayesdel import add_bayesdel_scores_to_built_file

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'test_files')


def test_add_bayesdel_scores(tmpdir):
    dest_path = os.path.join(tmpdir, 'output')

    runner = CliRunner()
    result = runner.invoke(add_bayesdel_scores_to_built_file.main,
                           [os.path.join(TEST_DATA_DIR, 'bayesdel', 'output.13.qc.vcf'),
                            os.path.join(TEST_DATA_DIR, 'bayesdel', 'output.17.qc.vcf'),
                            '--output', dest_path,
                            '--built-tsv', os.path.join(TEST_DATA_DIR, 'built_with_vr_ids_test_file.tsv')])

    assert result.exit_code == 0, result.stderr

    df = pd.read_csv(dest_path, sep='\t', na_values=['-'])
    assert 'BayesDel_nsfp33a_noAF' in df.columns
    np.testing.assert_equal(df['BayesDel_nsfp33a_noAF'].values, [-0.91, -0.95, np.nan, -0.12, np.nan])
