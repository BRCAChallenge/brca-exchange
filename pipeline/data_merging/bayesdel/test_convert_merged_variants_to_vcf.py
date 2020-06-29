import os

from data_merging.bayesdel import convert_merged_variants_to_vcf

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'test_files')

def test_convert_merged_to_vcf(tmp_path):
    vcf_path = os.path.join(TEST_DATA_DIR, 'built_with_vr_ids_test_file.tsv')

    output_path = tmp_path / 'file.vcf'
    convert_merged_variants_to_vcf.convert_merged_to_vcf(vcf_path, output_path)

    with open(output_path, 'r') as f:
        lines = f.readlines()

    # TOOD: more checks
    assert len(lines) == 9
