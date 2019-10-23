
import pytest
from common import variant_utils
from common.hgvs_utils import HgvsWrapper

variants = ['chr17:g.43094426:C>T', 'chr17:g.43071235:CCCT>C', 'chr13:g.32371231:T>TC',
            'chr13:g.32319239:CTCCAATAATATTCAAAGAGCAA>C', 'chr17:g.43125412:GCT>GAG']


@pytest.mark.parametrize("v", variants)
def test_variant_simple_roundtripping(v):
    assert str(variant_utils.VCFVariant.from_str(v)) == v


@pytest.mark.parametrize("v", variants)
def test_variant_hgvs_roundtripping(v):
    hgvs_wrap = HgvsWrapper()

    hgvs = variant_utils.VCFVariant.from_str(v).to_hgvs_obj(hgvs_wrap.contig_maps[hgvs_wrap.GRCh38_Assem])

    assert str(variant_utils.VCFVariant.from_hgvs_obj(hgvs)) == v
