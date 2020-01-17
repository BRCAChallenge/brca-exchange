import pandas as pd

from .utils import ChrInterval, build_interval_trees_by_chr

HGVS_CDNA_DEFAULT_AC = 'hgvs_cdna_default_ac'
SYMBOL_COL = 'symbol'
SYNONYM_AC_COL = 'synonyms_ac_col'


def load_config(path):
    '''
    Loads gene metadata from file into a dataframe

    :param path: config file path
    :return: dataframe
    '''
    df = pd.read_csv(path, sep=',', header=0, na_values='-')

    # allow for '-' in legacy variants to avoid duplicated information in the file,
    # in case the boundaries coincide
    df['start_hg38_legacy_variants'] = df['start_hg38_legacy_variants'].fillna(df['start_hg38']).astype(int)
    df['end_hg38_legacy_variants'] = df['end_hg38_legacy_variants'].fillna(df['end_hg38']).astype(int)

    return df.set_index(SYMBOL_COL, drop=False)


def get_genome_regions_symbol_dict(gene_config_df, start_col='start_hg38', end_col='end_hg38'):
    '''
    Build a dictionary from chromosome to a interval tree with the gene symbol
    as payload.

    Used to return the gene symbol (e.g. BRCA1) of a given position.

    :param gene_config_df:
    :return: dict[int, IntervalTree]
    '''

    def interval_tree_mapper(c, s, e):
        return gene_regions_dict[ChrInterval(c, s, e)]['symbol']

    gene_regions_dict = extract_gene_regions_dict(gene_config_df, start_col, end_col)

    return build_interval_trees_by_chr(gene_regions_dict.keys(),
                                    interval_tree_mapper)


def extract_gene_regions_dict(gene_config_df, start_col='start_hg38', end_col='end_hg38'):
    '''
    Prepares dict used in get_genome_regions_symbol_dict

    :param gene_config_df: gene metadata dataframe
    :return: dict[ChrInterval, dict[str, object]]
    '''

    # a[2]+1 -> end_hg38 position is inclusive, but in seq_utils end position is treated exclusive
    return {ChrInterval(a[0], a[1], a[2] + 1): {'symbol': a[3]} for a
            in
            gene_config_df.loc[:,
            ['chr', start_col, end_col, 'symbol']].values}
