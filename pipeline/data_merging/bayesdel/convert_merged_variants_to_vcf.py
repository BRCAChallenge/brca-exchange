import pandas as pd
import click

# TODO: use variant_utils?
def sp(s):
    return s.split(':')[0].lstrip('chr'), int(s.split(':')[1].lstrip('g.')), s, s.split(':')[2].split('>')[0], s.split(':')[2].split('>')[1]


def convert_merged_to_vcf(path, vcf_path):
    coord_col = 'Genomic_Coordinate_hg38'

    df = pd.read_csv(path, sep='\t')

    s = df.apply(lambda d: sp(d[coord_col]), axis=1)

    df_vcf = pd.DataFrame.from_dict(dict(zip(s.index, s.values))).T
    df_vcf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT']
    df_vcf['QUAL'] = df_vcf['FILTER'] = df_vcf['ID'] = '.'
    df_vcf['INFO'] = '.' #'Dummy1=nan;Dummy2=nan'

    df_vcf = df_vcf[~df_vcf['ALT'].apply(lambda s: any(x not in 'ACTG' for x in s))]

    with open(vcf_path, 'w') as f:
        # TODO split string up
        #header = '##fileformat=VCFv4.0\n##source=brcaexchange\n##reference=GRCh38\n##INFO=<ID=Dummy1,Number=.,Type=String,Description="">\n##INFO=<ID=Dummy2,Number=.,Type=String,Description="">\n'
        header = '##fileformat=VCFv4.0\n##source=brcaexchange\n##reference=GRCh38\n'
        f.write(header)
        f.write("#CHROM POS ID REF ALT QUAL FILTER INFO\n".replace(' ', '\t'))

    df_vcf.to_csv(vcf_path, sep='\t', mode='a', index=False, header=None)

# TODO use click or argparse?

@click.command()
@click.argument('merged_tsv', type=click.Path(readable=True))
@click.argument('output_vcf', type=click.Path(writable=True))
def main(merged_tsv, output_vcf):
    convert_merged_to_vcf(merged_tsv, output_vcf)


if __name__ == "__main__":
    main()
