import os
import pandas as pd

from data_merging.bayesdel import convert_merged_variants_to_vcf
from common import vcf_files_helper

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'test_files')


def test_convert_merged_to_vcf(tmp_path):
    vcf_path = os.path.join(TEST_DATA_DIR, 'built_with_vr_ids_test_file.tsv')

    output_path = tmp_path / 'file.vcf'
    convert_merged_variants_to_vcf.convert_merged_to_vcf(vcf_path, output_path)

    df = vcf_files_helper.read_vcf_as_dataframe(output_path)

    assert df[0].is_monotonic_increasing and df[1].is_monotonic_increasing

    expected_df = pd.read_csv(os.path.join(TEST_DATA_DIR, 'bayesdel', 'expected_merged_variants_vcf.csv'))
    expected_df.columns = [int(c) for c in expected_df.columns]

    pd.testing.assert_frame_equal(df, expected_df)
