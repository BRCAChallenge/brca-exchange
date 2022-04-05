import io
import os
import click
import pandas as pd
import pdb

coord_col = 'genomics_coord'

VCF_INFO_COL = 7
SPLICE_COLS = []

def read_vcf_as_dataframe_and_extract_column_names(path):
    with open(path, 'r') as f:
        lines = []
        for l in f:
            if l.startswith('##INFO'):
                info_cols = l.split('Format:')[1].split('|')
                for col in info_cols:
                    SPLICE_COLS.append(col.replace('">\n', '').strip())
            elif not l.startswith('##'):
                lines.append(l)
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


def spliceai_results_as_df(vcf):
    """
    Reading SpliceAI vcf and extracting useful information

    :return: dataframe with one column per variable of interest + join key
    """
    # spliceai_vcf_df = vcf_files_helper.read_vcf_as_dataframe(vcf).reset_index(drop=True)
    spliceai_vcf_df = read_vcf_as_dataframe_and_extract_column_names(vcf)

    # processing VCF INFO col, generate a dict variable_name -> value out of the info field for every record
    info_dict = spliceai_vcf_df.iloc[:, VCF_INFO_COL].str.split('|')
    # .apply(lambda l: {s.split('=')[0] : '='.join(s.split('=')[1:]) for s in l})

    # generating dataframe out of dict, one column per variable
    df_spliceai_props = pd.DataFrame.from_records(info_dict.values, index=info_dict.index)

    # update columns names
    df_spliceai_props.columns = SPLICE_COLS

    # join back to vcf dataframe
    df_ret = spliceai_vcf_df.merge(df_spliceai_props, how='inner', left_index=True, right_index=True)

    # calculate a coordinate representation to join with built_tsv
    df_ret[coord_col] = 'chr' + df_ret.iloc[:, 0].astype(str) + ":g." + df_ret.iloc[:, 1].astype(
        str) + ':' + df_ret.iloc[:, 3].astype(
        str) + ">" + df_ret.iloc[:, 4].astype(str)

    # dropping VCF columns, leaving columns for variable of interest + join field
    return df_ret.drop(columns=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'])


@click.command()
@click.option('--vcf', required=True, type=click.Path(readable=True))
@click.option('--built-tsv', required=True, type=click.Path(readable=True))
@click.option('--output', required=True, type=click.Path(writable=True))
def main(vcf, built_tsv, output):
    df_spliceai = spliceai_results_as_df(vcf)
    df = pd.read_csv(built_tsv, sep='\t', keep_default_na=False)

    df_merged = df.merge(df_spliceai, left_on='Genomic_Coordinate_hg38', right_on=coord_col, how='left')

    pdb.set_trace()

    # drop join key and write
    (df_merged.drop(columns=[coord_col]).
     to_csv(output, sep='\t', index=False, na_rep='-'))


if __name__ == "__main__":
    main()
