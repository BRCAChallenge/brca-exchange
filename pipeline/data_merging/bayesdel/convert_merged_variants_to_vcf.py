import click
import pandas as pd


def _coord_to_vcf_cols(coord):
    return (coord.split(':')[0].lstrip('chr'),  # chromosome col
            int(coord.split(':')[1].lstrip('g.')),  # pos col
            coord,  # id col
            coord.split(':')[2].split('>')[0],  # ref col
            coord.split(':')[2].split('>')[1])  # alt col


def convert_merged_to_vcf(path, vcf_path):
    coord_col = 'Genomic_Coordinate_hg38'

    df = pd.read_csv(path, sep='\t')

    vcf_cols = df.apply(lambda d: _coord_to_vcf_cols(d[coord_col]), axis=1)

    df_vcf = pd.DataFrame.from_dict(dict(zip(vcf_cols.index, vcf_cols.values))).T
    df_vcf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT']
    df_vcf['QUAL'] = df_vcf['FILTER'] = df_vcf['ID'] = '.'
    df_vcf['INFO'] = '.'

    # filtering out variants with unambiguous bases
    df_vcf = df_vcf[~df_vcf['ALT'].apply(lambda s: any(x not in 'ACTG' for x in s))]

    with open(vcf_path, 'w') as f:
        header = '##fileformat=VCFv4.0\n##source=brcaexchange\n##reference=GRCh38\n'
        f.write(header)
        f.write("#CHROM POS ID REF ALT QUAL FILTER INFO\n".replace(' ', '\t'))

    (df_vcf.sort_values(['CHROM', 'POS']).
            to_csv(vcf_path, sep='\t', mode='a', index=False, header=None))


@click.command()
@click.argument('merged_tsv', type=click.Path(readable=True))
@click.argument('output_vcf', type=click.Path(writable=True))
def main(merged_tsv, output_vcf):
    convert_merged_to_vcf(merged_tsv, output_vcf)


if __name__ == "__main__":
    main()
