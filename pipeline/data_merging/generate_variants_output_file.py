from pathlib import Path
from typing import Union, Dict

import click
import numpy as np
import pandas as pd

from common import hgvs_utils
from common import utils
from common import variant_utils

FIELD_NAME_FIELD = 'field_name'
CATEGORY_FIELD = 'category'
VARIANTS_OUTPUT_MAPPING_FIELD = 'variants_output_mapping'

NEW_COLUMN_NAME_FIELD = 'new_column_name'

TRANSFORM_TAG = "#TRANSFORM"
SKIP_TAG = "#SKIP"
COPY_TAG = "#COPY"

TRANSFORM_FIELDS = frozenset(['Alt', 'Pos', 'Ref'])

NORMALIZED_GENOMIC_COORD_FIELD = 'Genomic_HGVS_38'

GENOMIC_VCF38_FIELD = 'genomic_vcf38'


def _handle_transform_fields(df, df_metadata):
    # verifying metadata is in sync with code
    transform_fields = set(
        df_metadata[df_metadata[VARIANTS_OUTPUT_MAPPING_FIELD].str.upper() == TRANSFORM_TAG][FIELD_NAME_FIELD])
    assert TRANSFORM_FIELDS == transform_fields, \
        f"Expected {transform_fields} to be handled specially according to metadata file, but in code there is {TRANSFORM_FIELDS}."

    hgvs_wrapper = hgvs_utils.HgvsWrapper().get_instance()

    hgvs_obj = df[NORMALIZED_GENOMIC_COORD_FIELD].apply(
        lambda v: hgvs_wrapper.hgvs_parser.parse(v) if not isinstance(v, float) else v)

    # extracting alt, pos and ref from normalized hg38 coordinates in internal representation, e.g. chr13:g.32314514:C>T
    # in the original file, alt pos and ref are taken from the source data and may hence not be normalized correctly
    vcf_coords = [variant_utils.VCFVariant.from_hgvs_obj(h) if not isinstance(h, float) else h for h in hgvs_obj]

    # adding new field
    df.loc[:, GENOMIC_VCF38_FIELD] = pd.Series([str(v) if not isinstance(v, float) else v for v in vcf_coords])

    # transforming existing fields
    df.loc[:, 'Alt'] = pd.Series([v.alt if not isinstance(v, float) else v for v in vcf_coords])
    df.loc[:, 'Pos'] = pd.Series([v.pos if not isinstance(v, float) else v for v in vcf_coords])
    df.loc[:, 'Ref'] = pd.Series([v.ref if not isinstance(v, float) else v for v in vcf_coords])

    return df


def _determine_column_order(df_tmp, df_metadata):
    df_metadata_relevant_cols = df_metadata.loc[df_metadata[NEW_COLUMN_NAME_FIELD].isin(set(df_tmp.columns)),
                                                [CATEGORY_FIELD, NEW_COLUMN_NAME_FIELD]]

    # principally sorting by (category, column name)
    cols_sorted = list(df_metadata_relevant_cols.sort_values([CATEGORY_FIELD, NEW_COLUMN_NAME_FIELD])[
                           NEW_COLUMN_NAME_FIELD])

    # except putting a few key fields in front
    column_list_manual = ['gene_symbol', 'genomic_vcf38', 'pos', 'ref', 'alt', 'genomic_hgvs_38', 'cdna', 'protein', 'source',
                          'clinical_significance_enigma']

    column_list_final = column_list_manual + [c for c in cols_sorted if c not in set(column_list_manual)]

    assert set(df_tmp.columns) == set(column_list_final), f"Columns don't match, " \
                                                          f"missing in final columns: {set(df_tmp.columns) -  set(column_list_final)}, " \
                                                          f"additions in final columns: {set(column_list_final) -  set(df_tmp.columns)} "

    return column_list_final


def _write_field_metadata(df_metadata: pd.DataFrame, output_metadata: Union[Path, str]):
    df_metadata_var_output = (df_metadata[df_metadata[VARIANTS_OUTPUT_MAPPING_FIELD] != SKIP_TAG]
                              .drop(columns=[VARIANTS_OUTPUT_MAPPING_FIELD, FIELD_NAME_FIELD])
                              .rename(columns={NEW_COLUMN_NAME_FIELD: FIELD_NAME_FIELD}))

    columns_ordered = [FIELD_NAME_FIELD] + [c for c in df_metadata_var_output.columns if c != FIELD_NAME_FIELD]
    df_metadata_var_output = df_metadata_var_output[columns_ordered].sort_values(FIELD_NAME_FIELD)

    utils.write_dataframe_as_tsv(df_metadata_var_output, output_metadata)


def _get_renaming_map(df_metadata: pd.DataFrame) -> Dict[str, str]:
    is_skip = df_metadata[VARIANTS_OUTPUT_MAPPING_FIELD] == SKIP_TAG
    is_copy_or_transform = (df_metadata[VARIANTS_OUTPUT_MAPPING_FIELD] == COPY_TAG) | (
            df_metadata[VARIANTS_OUTPUT_MAPPING_FIELD] == TRANSFORM_TAG)
    is_rename = (~is_skip) & (~is_copy_or_transform)

    df_metadata.loc[is_skip, NEW_COLUMN_NAME_FIELD] = ""
    df_metadata.loc[is_copy_or_transform, NEW_COLUMN_NAME_FIELD] = df_metadata.loc[
        is_copy_or_transform, FIELD_NAME_FIELD].str.lower()
    df_metadata.loc[is_rename, NEW_COLUMN_NAME_FIELD] = df_metadata.loc[
        is_rename, VARIANTS_OUTPUT_MAPPING_FIELD].str.lower()

    return {old_name: new_name for (_, (old_name, new_name)) in
            df_metadata[[FIELD_NAME_FIELD, NEW_COLUMN_NAME_FIELD]].iterrows() if new_name}


@click.command()
@click.argument('build_file', type=click.Path(readable=True))
@click.argument('field_metadata', type=click.Path(readable=True))
@click.argument('field_metadata_additional', type=click.Path(readable=True))
@click.argument('output', type=click.Path(writable=True))
@click.argument('output_metadata', type=click.Path(writable=True))
def main(build_file, field_metadata, field_metadata_additional, output, output_metadata):
    df = utils.read_tsv_as_dataframe(build_file)

    df_metadata = utils.read_tsv_as_dataframe(field_metadata)

    # verifying input dataframe and metadata are in sync
    assert set(df.columns) == set(df_metadata[FIELD_NAME_FIELD]), f"Field in metadata and columns in input data don't match. " \
                                                                  f"Missing in metadata: {set(df.columns) - set(df_metadata[FIELD_NAME_FIELD])} " \
                                                                  f"Missing in dataframe: {set(df_metadata[FIELD_NAME_FIELD]) - set(df.columns)}."

    # some fields are not copied from the original data as is, but are transformed
    df_transformed = _handle_transform_fields(df, df_metadata)

    df_metadata_all = pd.concat([df_metadata, utils.read_tsv_as_dataframe(field_metadata_additional)])
    # rename columns
    col_map_renaming = _get_renaming_map(df_metadata_all)
    column_list = col_map_renaming.values()

    df_renamed_cols = df_transformed.rename(columns=col_map_renaming).loc[:, column_list]

    column_list_final = _determine_column_order(df_renamed_cols, df_metadata_all)

    # project and reorder columns
    df_out = df_renamed_cols.loc[:, column_list_final]

    df_out.replace('', np.nan, inplace=True)

    assert df.shape[0] == df_out.shape[0], "Number of rows differs!"

    utils.write_dataframe_as_tsv(df_out, output)

    # write out field metadata file tailored to the variants output file
    _write_field_metadata(df_metadata_all, output_metadata)


if __name__ == "__main__":
    main()
