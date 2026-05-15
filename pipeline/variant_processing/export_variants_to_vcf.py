"""
Export all GRCh38 variants from the pipeline DB to a VCF file.

Produces a minimal VCF (no INFO) suitable as input to SpliceAI.
Only variants whose REF and ALT consist entirely of ACGT bases are included.
"""

import click
import psycopg2
import psycopg2.extras

VCF_HEADER = """\
##fileformat=VCFv4.0
##source=brcaexchange
##reference=GRCh38
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=8,length=146364022>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=18,length=78077248>
##contig=<ID=19,length=59128983>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=22,length=51304566>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

_ACGT = set('ACGT')


def _is_clean(seq):
    return bool(seq) and all(c in _ACGT for c in seq.upper())


@click.command()
@click.option('--db-url', default='postgresql://postgres:postgres@localhost/storage.pg',
              envvar='PIPELINE_DB_URL', show_default=True)
@click.option('--schema', default='pipeline', show_default=True)
@click.option('--output', required=True, type=click.Path(writable=True),
              help='Output VCF path')
def main(db_url, schema, output):
    conn = psycopg2.connect(db_url, options=f'-c search_path={schema}')
    try:
        with conn.cursor() as cur:
            cur.execute("""
                SELECT chr, pos, ref, alt
                FROM genomic_coordinates
                WHERE assembly = 'GRCh38'
                ORDER BY chr, pos::bigint
            """)
            rows = cur.fetchall()
    finally:
        conn.close()

    written = 0
    with open(output, 'w') as fh:
        fh.write(VCF_HEADER)
        for chr_, pos, ref, alt in rows:
            if not (_is_clean(ref) and _is_clean(alt)):
                continue
            chr_bare = chr_[3:] if chr_.startswith('chr') else chr_
            fh.write(f'{chr_bare}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n')
            written += 1

    print(f'Exported {written} variants to {output}')


if __name__ == '__main__':
    main()
