
import pytest
from common import seq_utils, variant_utils
from common.hgvs_utils import HgvsWrapper

variants = ['chr17:g.43094426:C>T',
            'chr17:g.43071235:CCCT>C',
            'chr13:g.32371231:T>TC',
            'chr13:g.32319239:CTCCAATAATATTCAAAGAGCAA>C',
            'chr17:g.43125412:GCT>GAG']


@pytest.fixture(scope="module")
def wrapper():
    return HgvsWrapper()


@pytest.fixture(scope="module")
def seq_fetcher():
    return seq_utils.SeqRepoWrapper.get_instance()


@pytest.mark.parametrize("v", variants)
def test_variant_simple_roundtripping(v):
     assert str(variant_utils.VCFVariant.from_str(v)) == v


@pytest.mark.parametrize("v", variants)
def test_variant_hgvs_roundtripping(v):
    hgvs_wrap = HgvsWrapper()
    hgvs = variant_utils.VCFVariant.from_str(v).to_hgvs_obj(hgvs_wrap.contig_maps[hgvs_wrap.GRCh38_Assem])
    assert str(variant_utils.VCFVariant.from_hgvs_obj(hgvs)) == v


@pytest.mark.parametrize("v,expected",
                         [('NC_000017.11:g.43091923C>T', 'chr17:g.43091923:C>T'),
                          ('NC_000017.11:g.43086110_43086115dupACACAC', 'chr17:g.43086109:T>TACACAC'),
                          ('NC_000017.11:g.43098661_43098662insT',  'chr17:g.43098661:C>CT'),
                          ('NC_000017.11:g.43091784_43091817del', 'chr17:g.43091783:CGGTAGCAACGGTGCTATGCCTAGTAGACTGAGAA>C')
                        ])
def test_hgvs_to_vcf(wrapper, v, expected, seq_fetcher):
    hgvs_var = wrapper.hgvs_parser.parse(v)
    v = variant_utils.VCFVariant.from_hgvs_obj(hgvs_var, seq_fetcher)

    assert str(v) == expected
