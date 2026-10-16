import os

import pytest

from insilico import make_spliceai_annotation as msa

ANNOTATION = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'spliceai_annotations', 'brca_mane_grch38.txt')

# Two exons on each strand, in GTF order (BRCA2 ascending, BRCA1 descending),
# to check that we sort into ascending genomic order regardless.
GTF = '''\
#!annotation-source MANE test fixture
chr13\tx\ttranscript\t100\t400\t.\t+\t.\tgene_name "BRCA2"; transcript_id "ENST2.1"; tag "MANE_Select"; db_xref "RefSeq:NM_2.4";
chr13\tx\texon\t100\t150\t.\t+\t.\tgene_name "BRCA2"; transcript_id "ENST2.1"; tag "MANE_Select";
chr13\tx\texon\t300\t400\t.\t+\t.\tgene_name "BRCA2"; transcript_id "ENST2.1"; tag "MANE_Select";
chr17\tx\ttranscript\t900\t2000\t.\t-\t.\tgene_name "BRCA1"; transcript_id "ENST1.9"; tag "MANE_Select"; db_xref "RefSeq:NM_1.4";
chr17\tx\texon\t1800\t2000\t.\t-\t.\tgene_name "BRCA1"; transcript_id "ENST1.9"; tag "MANE_Select";
chr17\tx\texon\t900\t1000\t.\t-\t.\tgene_name "BRCA1"; transcript_id "ENST1.9"; tag "MANE_Select";
chr17\tx\ttranscript\t900\t2000\t.\t-\t.\tgene_name "BRCA1"; transcript_id "ENST1.8"; tag "MANE_Plus_Clinical";
chr7\tx\ttranscript\t10\t20\t.\t+\t.\tgene_name "OTHER"; transcript_id "ENST9.1"; tag "MANE_Select";
'''


@pytest.fixture
def gtf(tmp_path):
    path = tmp_path / 'mane.gtf'
    path.write_text(GTF)
    return str(path)


def read(path):
    with open(path) as fp:
        header = fp.readline().rstrip('\n').split('\t')
        return header, [dict(zip(header, line.rstrip('\n').split('\t')))
                        for line in fp]


def test_writes_spliceai_columns_with_zero_based_starts(gtf, tmp_path):
    out = str(tmp_path / 'ann.txt')
    msa.write_annotation(msa.read_gtf(gtf, ['BRCA1', 'BRCA2'], 'MANE_Select', False),
                         ['BRCA1', 'BRCA2'], out)
    header, rows = read(out)
    assert header == ['#NAME', 'CHROM', 'STRAND', 'TX_START', 'TX_END',
                      'EXON_START', 'EXON_END']
    brca1, brca2 = rows
    # GTF is 1-based inclusive; starts drop by one, ends are unchanged.
    assert (brca1['CHROM'], brca1['STRAND']) == ('17', '-')
    assert (brca1['TX_START'], brca1['TX_END']) == ('899', '2000')
    # Ascending genomic order on the minus strand too, comma-terminated.
    assert brca1['EXON_START'] == '899,1799,'
    assert brca1['EXON_END'] == '1000,2000,'
    assert (brca2['TX_START'], brca2['TX_END']) == ('99', '400')
    assert brca2['EXON_START'] == '99,299,'


def test_keeps_chr_prefix_only_when_asked(gtf, tmp_path):
    out = str(tmp_path / 'ann.txt')
    msa.write_annotation(msa.read_gtf(gtf, ['BRCA1'], 'MANE_Select', True),
                         ['BRCA1'], out)
    assert read(out)[1][0]['CHROM'] == 'chr17'


def test_ignores_other_genes_and_other_tags(gtf):
    records = msa.read_gtf(gtf, ['BRCA1', 'BRCA2'], 'MANE_Select', False)
    assert set(records) == {'BRCA1', 'BRCA2'}
    assert records['BRCA1']['transcript'] == 'ENST1.9'


def test_missing_gene_is_an_error(gtf):
    with pytest.raises(SystemExit, match='NOTAGENE'):
        msa.read_gtf(gtf, ['NOTAGENE'], 'MANE_Select', False)


def test_checked_in_table_matches_the_mane_select_transcripts():
    """The table SpliceAI is actually run with, guarded against hand edits."""
    header, rows = read(ANNOTATION)
    by_gene = {row['#NAME']: row for row in rows}
    assert set(by_gene) == {'BRCA1', 'BRCA2'}
    # NM_007294.4 / ENST00000357654.9 and NM_000059.4 / ENST00000380152.8.
    assert (by_gene['BRCA1']['CHROM'], by_gene['BRCA1']['STRAND']) == ('17', '-')
    assert (by_gene['BRCA1']['TX_START'], by_gene['BRCA1']['TX_END']) == \
        ('43044294', '43125364')
    assert (by_gene['BRCA2']['CHROM'], by_gene['BRCA2']['STRAND']) == ('13', '+')
    assert (by_gene['BRCA2']['TX_START'], by_gene['BRCA2']['TX_END']) == \
        ('32315507', '32400268')
    assert len(by_gene['BRCA1']['EXON_START'].split(',')) - 1 == 23
    assert len(by_gene['BRCA2']['EXON_START'].split(',')) - 1 == 27


def test_checked_in_table_parses_the_way_spliceai_parses_it():
    annotator = pytest.importorskip('spliceai.utils').Annotator
    import pandas

    frame = pandas.read_csv(ANNOTATION, sep='\t', dtype={'CHROM': object})
    assert list(frame['#NAME']) == ['BRCA1', 'BRCA2']
    assert annotator  # the columns above are the ones Annotator.__init__ reads
