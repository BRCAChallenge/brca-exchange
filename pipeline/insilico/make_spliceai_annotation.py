#!/usr/bin/env python
"""make_spliceai_annotation: build a SpliceAI gene annotation table for the
reference transcripts of the genes we score.

SpliceAI ships an annotation table (``-A grch38``) taken from GENCODE V24
canonical, and that table is not just used to name genes and mask scores: in
``spliceai.utils.get_delta_scores`` every base of the input window that falls
outside the annotated transcript is overwritten with N before one-hot
encoding.  Since the model's receptive field is 5,000 bases either side of
each scored position, a transcript whose ends are in the wrong place changes
the scores of every variant within ~5 kb of them, and variants outside the
span get no score at all.

The packaged BRCA1 record spans 43,045,629-43,125,483, which is 1,334 bases
short of the MANE Select transcript NM_007294.4 / ENST00000357654.9
(43,044,295-43,125,364) at the 3' end.  This table is built from MANE instead,
so that we score against the same reference transcripts the ENIGMA BRCA1/2
VCEP classifies against.

The output is SpliceAI's own format: a tab-separated table with the columns
NAME, CHROM, STRAND, TX_START, TX_END, EXON_START, EXON_END, where the start
coordinates are 0-based and the exon lists are comma-terminated.  Note that
SpliceAI reads it with a bare ``pandas.read_csv`` and has no notion of comment
lines, so provenance lives in Readme.md next to the generated file rather than
in the file itself.

Usage:

    python make_spliceai_annotation.py \
        -g /path/to/MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz \
        -o insilico/spliceai_annotations/brca_mane_grch38.txt
"""

import argparse
import gzip
import re

DEFAULT_GENES = ['BRCA1', 'BRCA2']
DEFAULT_TAG = 'MANE_Select'

_ATTR_RE = re.compile(r'(\S+) "([^"]*)"')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--mane_gtf", required=True,
                        help="MANE genomic GTF (plain or gzipped)")
    parser.add_argument("-o", "--output", required=True,
                        help="Pathname for the SpliceAI annotation table")
    parser.add_argument("-n", "--genes", default=",".join(DEFAULT_GENES),
                        help="Comma-separated gene symbols to include")
    parser.add_argument("-t", "--tag", default=DEFAULT_TAG,
                        help="GTF tag identifying the reference transcript")
    parser.add_argument("--keep_chr_prefix", action="store_true",
                        help="Emit chr17 rather than 17 in the CHROM column")
    return parser.parse_args()


def attributes(field):
    return dict(_ATTR_RE.findall(field))


def read_gtf(mane_gtf, genes, tag, keep_chr_prefix):
    """Return {gene: record} for the tagged transcript of each wanted gene."""
    wanted = set(genes)
    records = {}
    opener = gzip.open if mane_gtf.endswith('.gz') else open
    with opener(mane_gtf, 'rt') as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9 or cols[2] not in ('transcript', 'exon'):
                continue
            attrs = attributes(cols[8])
            gene = attrs.get('gene_name')
            if gene not in wanted or tag not in cols[8]:
                continue
            chrom = cols[0]
            if not keep_chr_prefix and chrom.startswith('chr'):
                chrom = chrom[3:]
            # GTF is 1-based inclusive; SpliceAI wants 0-based starts.
            start, end = int(cols[3]) - 1, int(cols[4])
            record = records.setdefault(gene, {
                'chrom': chrom,
                'strand': cols[6],
                'transcript': attrs.get('transcript_id'),
                'refseq': attrs.get('db_xref', ''),
                'exons': [],
            })
            if cols[2] == 'transcript':
                record['tx_start'] = start
                record['tx_end'] = end
            else:
                record['exons'].append((start, end))

    missing = wanted - set(records)
    if missing:
        raise SystemExit('No {} transcript found in {} for: {}'.format(
            tag, mane_gtf, ', '.join(sorted(missing))))
    return records


def write_annotation(records, genes, output):
    with open(output, 'w') as fp:
        fp.write('#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n')
        for gene in genes:
            record = records[gene]
            # SpliceAI indexes the exon lists positionally and takes a min()
            # over their union, so the order only has to be consistent; the
            # packaged table lists them in ascending genomic order on both
            # strands, and we follow it.
            exons = sorted(record['exons'])
            if exons[0][0] != record['tx_start'] or exons[-1][1] != record['tx_end']:
                raise SystemExit(
                    '{} exons {}-{} do not span the transcript {}-{}'.format(
                        gene, exons[0][0], exons[-1][1],
                        record['tx_start'], record['tx_end']))
            fp.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                gene, record['chrom'], record['strand'],
                record['tx_start'], record['tx_end'],
                ','.join(str(s) for s, _ in exons) + ',',
                ','.join(str(e) for _, e in exons) + ','))
            print('{}: {} {}:{}-{} ({} exons, {})'.format(
                gene, record['transcript'], record['chrom'],
                record['tx_start'] + 1, record['tx_end'],
                len(exons), record['refseq'] or 'no RefSeq xref'))


def main():
    args = parse_args()
    genes = [g.strip() for g in args.genes.split(',') if g.strip()]
    records = read_gtf(args.mane_gtf, genes, args.tag, args.keep_chr_prefix)
    write_annotation(records, genes, args.output)
    print('Wrote', args.output)


if __name__ == '__main__':
    main()
