"""
A few integration tests
"""

import pytest
from common import hgvs_utils
from common import variant_utils


@pytest.fixture()
def wrapper():
    return hgvs_utils.HgvsWrapper()


@pytest.fixture()
def hgvs_variant(wrapper):
    return (variant_utils.VCFVariant.from_str('chr13:g.32316477:AAG>A').
            to_hgvs_obj(wrapper.contig_maps[hgvs_utils.HgvsWrapper.GRCh38_Assem]))


def test_to_cdna(wrapper, hgvs_variant):
    cdna = wrapper.to_cdna(hgvs_variant)
    assert str(cdna) == 'NM_000059.3:c.17_19delinsA' # normalized 'NM_000059.3:c.22_23del'


def test_to_protein(wrapper, hgvs_variant):
    protein = wrapper.to_protein(wrapper.to_cdna(hgvs_variant))
    assert str(protein) == 'NP_000050.2:p.(Arg8AlafsTer5)'


def test_normalize(wrapper, hgvs_variant):
    obj = wrapper.normalizing(hgvs_variant)
    assert str(obj) == 'NC_000013.11:g.32316482_32316483del'
