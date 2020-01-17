import clinvar

from common import hgvs_utils, variant_utils
import xml.etree.ElementTree as ET


def test_simple_genomic_coordinate_extraction():
    sample_measure = """
      <Measure Type="Deletion" ID="24358">
        <Name>
          <ElementValue Type="Preferred">NM_000059.3(BRCA2):c.6591_6592del (p.Glu2198fs)</ElementValue>
        </Name>
        <SequenceLocation Assembly="GRCh38" AssemblyAccessionVersion="GCF_000001405.38" AssemblyStatus="current" Chr="13" Accession="NC_000013.11" start="32340946" stop="32340947" display_start="32340946" display_stop="32340947" variantLength="2" positionVCF="32340945" referenceAlleleVCF="CTG" alternateAlleleVCF="C"/>
        <SequenceLocation Assembly="GRCh37" AssemblyAccessionVersion="GCF_000001405.25" AssemblyStatus="previous" Chr="13" Accession="NC_000013.10" start="32915083" stop="32915084" display_start="32915083" display_stop="32915084" variantLength="2" positionVCF="32915082" referenceAlleleVCF="CTG" alternateAlleleVCF="C"/>
     </Measure>
    """

    measure_el = ET.fromstring(sample_measure)

    genomic_coords = clinvar.extract_genomic_coordinates_from_measure(measure_el)

    assert genomic_coords[hgvs_utils.HgvsWrapper.GRCh38_Assem] == variant_utils.VCFVariant(13, 32340945, "CTG", "C")


def test_genomic_coordinate_extraction_from_NM():
    sample_measure = """
    <Measure Type="Deletion" ID="24358">
        <Name>
          <ElementValue Type="Preferred">NM_000059.3(BRCA2):c.6591_6592del (p.Glu2198fs)</ElementValue>
        </Name>
    </Measure>
    """
    measure_el = ET.fromstring(sample_measure)

    genomic_coords = clinvar.extract_genomic_coordinates_from_measure(measure_el)

    assert genomic_coords[hgvs_utils.HgvsWrapper.GRCh38_Assem] == variant_utils.VCFVariant(13, 32340945, "CTG", "C")


def test_preprocess_element_value():
    assert clinvar._preprocess_element_value('NM_000059.3(BRCA2):c.6591_6592del (p.Glu2198fs)') == 'NM_000059.3(BRCA2):c.6591_6592del'







