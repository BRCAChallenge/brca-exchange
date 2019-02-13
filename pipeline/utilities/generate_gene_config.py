import click
import mygene
import pandas as pd

mg = mygene.MyGeneInfo()


def _process_response(r):
    return {
        'symbol': r['symbol'],
        'entrez_id': r['entrezgene'],
        'ensembl_id': r['ensembl']['gene'],
        'chr': r['genomic_pos']['chr'],
        'start_hg37': r['genomic_pos_hg19']['start'],
        'end_hg37': r['genomic_pos_hg19']['end'],
        'start_hg38': r['genomic_pos']['start'],
        'end_hg38': r['genomic_pos']['end']
    }


def _get_entrez_ids(gene_symbols):
    entrez_ids = []
    for gene in gene_symbols:
        gene = gene.upper()

        ret = mg.query('symbol:{}'.format(gene), species='human')

        assert (ret['total'] == 1)

        entrez_ids.append(ret['hits'][0]["entrezgene"])

    return entrez_ids


def _generate_gene_df(gene_symbols):
    entrez_ids = _get_entrez_ids(gene_symbols)

    ret = mg.getgenes(entrez_ids,
                      fields="symbol,genomic_pos_hg19,genomic_pos,ensembl,entrezgene")

    df = pd.DataFrame.from_records(
        [_process_response(r) for r in ret]).sort_values('symbol')

    df.sort_values('symbol')

    # putting symbol in front
    cols = ['symbol', 'entrez_id']
    cols.extend([c for c in df.columns if c not in cols])
    df = df[cols]

    return df


@click.command()
@click.argument('gene_symbols', nargs=-1)
@click.option("-o", help="output path")
def main(gene_symbols, o):
    df = _generate_gene_df(gene_symbols)

    df.to_csv(o, header=True, index=False)


if __name__ == '__main__':
    main()
