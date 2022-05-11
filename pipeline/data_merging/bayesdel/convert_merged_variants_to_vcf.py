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

    # only taking variants with unambiguous bases
    df_vcf = df_vcf[df_vcf['ALT'].apply(lambda s: all(x in 'ACTG' for x in s))]

    with open(vcf_path, 'w') as f:
        header = '##fileformat=VCFv4.0\n##source=brcaexchange\n##reference=GRCh38\n'
        f.write(header)
        f.write(("##contig=<ID=1,length=249250621>"
                 "\n##contig=<ID=2,length=243199373>"
                 "\n##contig=<ID=3,length=198022430>"
                 "\n##contig=<ID=4,length=191154276>"
                 "\n##contig=<ID=5,length=180915260>"
                 "\n##contig=<ID=6,length=171115067>"
                 "\n##contig=<ID=7,length=159138663>"
                 "\n##contig=<ID=8,length=146364022>"
                 "\n##contig=<ID=9,length=141213431>"
                 "\n##contig=<ID=10,length=135534747>"
                 "\n##contig=<ID=11,length=135006516>"
                 "\n##contig=<ID=12,length=133851895>"
                 "\n##contig=<ID=13,length=115169878>"
                 "\n##contig=<ID=14,length=107349540>"
                 "\n##contig=<ID=15,length=102531392>"
                 "\n##contig=<ID=16,length=90354753>"
                 "\n##contig=<ID=17,length=81195210>"
                 "\n##contig=<ID=18,length=78077248>"
                 "\n##contig=<ID=19,length=59128983>"
                 "\n##contig=<ID=20,length=63025520>"
                 "\n##contig=<ID=21,length=48129895>"
                 "\n##contig=<ID=22,length=51304566>\n"))
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
