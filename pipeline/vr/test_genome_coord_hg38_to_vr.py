import pytest
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.parser
import genomic_coord_hg38_to_vr as hg382vr
import re
import VR

@pytest.fixture
def hg38_mapper():
    """Returns a mapper to map strings to hg38 variant mapper"""
    return VR.Mapper()

def test_hg38_mapping_basic_variant(hg38_mapper):
    hg38_variant = hg38_mapper.map_variant("NM_000059.3:c.-200C>T")
    assert(str(hg38_variant) == "NC_000013.11:g.32315507C>T")


def test_hg38_mapping_out_of_bounds_variant(hg38_mapper):
    hg38_variant = hg38_mapper.map_variant("NM_000059.3:c.-1193C>T")
    assert(hg38_variant is None)


def test_hg38_mapping_deletion_variant(hg38_mapper):
    hg38_variant = hg38_mapper.map_variant("NM_000059.3:c.-59_-57delGAA")

    # use a regular expression match instead of string comparison, so that the
    # test doesn't depend on whether or not the deleted residues are part of
    # the name.
    match = re.search("^NC_000013.11:g.32315648_32315650del",
                      str(hg38_variant))
    assert(match is not None)


def test_hg38_mapping_insertion_variant(hg38_mapper):
    hg38_variant = hg38_mapper.map_variant("NM_000059.3:c.37_38insT")
    assert(str(hg38_variant) == "NC_000013.11:g.32316497_32316498insT")


@pytest.fixture
def translator():
    """Returns a translator to translate data to VR"""
    return VR.VRUtils("http://localhost:5000/seqrepo")


def test_VR_missense_variant(translator, hg38_mapper):
    """Test translation of a missense variant to VR"""
    vr_translation = translator.to_vr("NM_000059.3:c.-196G>A")
    vr = VR.VR(vr_translation)
    assert(vr.start() == 32315510 and vr.end() == 32315511
           and vr.alt() == "A")

def test_VR_deletion_variant(translator, hg38_mapper):
    """Test translation of a deletion variant to VR"""
    vr_translation = translator.to_vr("NM_000059.3:c.37_44delGAAATTTT")
    vr = VR.VR(vr_translation)
    assert(vr.start() == 32316496 and vr.end() == 32316504
           and vr.alt() == "")


def test_VR_indel_variant(translator, hg38_mapper):
    """Test translation of a variant with both an insertion and a deletion"""
    vr_translation = translator.to_vr("NM_000059.3:c.32_33delTTinsA")
    vr = VR.VR(vr_translation)
    assert(vr.start() == 32316491 and vr.end() == 32316493
           and vr.alt() == "A")
    
    


