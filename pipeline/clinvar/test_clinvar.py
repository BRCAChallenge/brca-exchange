from . import clinvar_common

from common import hgvs_utils, variant_utils
import xml.etree.ElementTree as ET


def test_simple_genomic_coordinate_extraction():
    sample_location = """
      <Location>
        <CytogeneticLocation>13q13.1</CytogeneticLocation>
        <SequenceLocation Assembly="GRCh38" AssemblyAccessionVersion="GCF_000001405.38" forDisplay="true" AssemblyStatus="cur\
rent" Chr="13" Accession="NC_000013.11" start="32345245" stop="32345245" display_start="32345245" display_stop="32345245" var\
iantLength="1" positionVCF="32345245" referenceAlleleVCF="A" alternateAlleleVCF="G"/>
        <SequenceLocation Assembly="GRCh37" AssemblyAccessionVersion="GCF_000001405.25" AssemblyStatus="previous" Chr="13" Ac\
cession="NC_000013.10" start="32919382" stop="32919382" display_start="32919382" display_stop="32919382" variantLength="1" po\
sitionVCF="32919382" referenceAlleleVCF="A" alternateAlleleVCF="G"/>
      </Location>
    """

    location_el = ET.fromstring(sample_location)

    genomic_coords = clinvar_common.extract_genomic_coordinates_from_location(location_el)

    assert genomic_coords[hgvs_utils.HgvsWrapper.GRCh38_Assem] == variant_utils.VCFVariant(13, 32345245, "A", "G")


def test_genomic_coordinate_extraction_from_NM():

    sample_name="NM_007294.4:c.6591_6592del"
    genomic_coords = clinvar_common.accession_to_genomic_coordinates(sample_name)

    assert genomic_coords[hgvs_utils.HgvsWrapper.GRCh38_Assem] == variant_utils.VCFVariant(17, 43044677, "AAT", "A")


def test_preprocess_element_value():
    assert clinvar_common._preprocess_element_value('NM_000059.3(BRCA2):c.6591_6592del (p.Glu2198fs)') == 'NM_000059.3(BRCA2):c.6591_6592del'







