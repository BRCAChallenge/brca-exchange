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

    # 3'UTR position (NM_007294.3's CDS is 5592bp, so a plain c.6591 position,
    # as opposed to c.*999, is out of bounds and would raise an HGVS error).
    sample_name="NM_007294.3:c.*999_*1000del"
    genomic_coords = clinvar_common.accession_to_genomic_coordinates(sample_name)

    assert genomic_coords[hgvs_utils.HgvsWrapper.GRCh38_Assem] == variant_utils.VCFVariant(17, 43044677, "AAT", "A")


def test_preprocess_element_value():
    assert clinvar_common._preprocess_element_value('NM_000059.3(BRCA2):c.6591_6592del (p.Glu2198fs)') == 'NM_000059.3(BRCA2):c.6591_6592del'


def test_is_bic_designation():
    assert clinvar_common.is_bic_designation('1294del41')
    assert clinvar_common.is_bic_designation('5277A>G')
    assert clinvar_common.is_bic_designation('999insA')
    assert not clinvar_common.is_bic_designation('p.(Leu392GlnfsTer6)')
    assert clinvar_common.is_bic_designation('S76X')
    assert clinvar_common.is_bic_designation('R245X')
    assert clinvar_common.is_bic_designation('E1703V')
    assert clinvar_common.is_bic_designation('S309=')
    assert not clinvar_common.is_bic_designation('NM_000059.4(BRCA2):c.67G>A')
    assert not clinvar_common.is_bic_designation('NM_007294.4:c.3541del')
    assert not clinvar_common.is_bic_designation('NC_000013.10:g.32900419G>T')
    assert not clinvar_common.is_bic_designation('LRG_292t1:c.2066_2069delGTAA')


def test_gene_chromosome_mismatch_rejected():
    # GeneList tags this variant as BRCA1, but its own coordinates are on
    # chromosome 13 (BRCA2's chromosome) -- a real ClinVar data-quality issue.
    sample_simple_allele = """
      <SimpleAllele AlleleID="12345">
        <Location>
          <SequenceLocation Assembly="GRCh38" Chr="13" positionVCF="32345245" referenceAlleleVCF="A" alternateAlleleVCF="G"/>
        </Location>
        <GeneList>
          <Gene Symbol="BRCA1"/>
        </GeneList>
      </SimpleAllele>
    """

    element = ET.fromstring(sample_simple_allele)
    v = clinvar_common.variant(element, gene_chromosomes={'BRCA1': 17, 'BRCA2': 13}, debug=False)

    assert v.geneSymbol is None


def test_gene_chromosome_match_accepted():
    sample_simple_allele = """
      <SimpleAllele AlleleID="12345">
        <Location>
          <SequenceLocation Assembly="GRCh38" Chr="17" positionVCF="43123988" referenceAlleleVCF="C" alternateAlleleVCF="T"/>
        </Location>
        <GeneList>
          <Gene Symbol="BRCA1"/>
        </GeneList>
      </SimpleAllele>
    """

    element = ET.fromstring(sample_simple_allele)
    v = clinvar_common.variant(element, gene_chromosomes={'BRCA1': 17, 'BRCA2': 13}, debug=False)

    assert v.geneSymbol == 'BRCA1'


def test_gene_chromosome_check_skipped_when_not_provided():
    # Without gene_chromosomes, behavior is unchanged: no cross-check happens.
    sample_simple_allele = """
      <SimpleAllele AlleleID="12345">
        <Location>
          <SequenceLocation Assembly="GRCh38" Chr="13" positionVCF="32345245" referenceAlleleVCF="A" alternateAlleleVCF="G"/>
        </Location>
        <GeneList>
          <Gene Symbol="BRCA1"/>
        </GeneList>
      </SimpleAllele>
    """

    element = ET.fromstring(sample_simple_allele)
    v = clinvar_common.variant(element, debug=False)

    assert v.geneSymbol == 'BRCA1'


def test_variant_bic_nomenclature_from_other_name_list():
    # OtherNameList mixes a legacy BIC-style designator with an HGVS-style
    # protein change synonym; only the former should be picked out.
    sample_simple_allele = """
      <SimpleAllele AlleleID="12345">
        <Location>
          <SequenceLocation Assembly="GRCh38" Chr="13" positionVCF="32345245" referenceAlleleVCF="A" alternateAlleleVCF="G"/>
        </Location>
        <OtherNameList>
          <Name>1294del41</Name>
          <Name>p.(Leu392GlnfsTer6)</Name>
        </OtherNameList>
      </SimpleAllele>
    """

    element = ET.fromstring(sample_simple_allele)
    v = clinvar_common.variant(element, debug=False)

    assert v.bic_nomenclature == '1294del41'
    assert 'p.(Leu392GlnfsTer6)' in v.synonyms


def test_variant_bic_nomenclature_absent():
    sample_simple_allele = """
      <SimpleAllele AlleleID="12345">
        <Location>
          <SequenceLocation Assembly="GRCh38" Chr="13" positionVCF="32345245" referenceAlleleVCF="A" alternateAlleleVCF="G"/>
        </Location>
        <OtherNameList>
          <Name>p.(Leu392GlnfsTer6)</Name>
        </OtherNameList>
      </SimpleAllele>
    """

    element = ET.fromstring(sample_simple_allele)
    v = clinvar_common.variant(element, debug=False)

    assert v.bic_nomenclature is None


# A BRCA1 VariationArchive whose two RCVs disagree with each other; the
# aggregate GermlineClassification is substituted in per test.
VARIATION_ARCHIVE_TEMPLATE = """
<VariationArchive VariationID="37394" VariationName="NM_007294.4(BRCA1):c.1175_1215del (p.Leu392fs)" Accession="VCV000037394" Version="12" DateLastUpdated="2025-03-01">
  <RecordStatus>current</RecordStatus>
  <ClassifiedRecord>
    <SimpleAllele AlleleID="46183" VariationID="37394">
      <GeneList><Gene Symbol="BRCA1"/></GeneList>
      <Location>
        <SequenceLocation Assembly="GRCh38" Chr="17" positionVCF="43094315" referenceAlleleVCF="TTG" alternateAlleleVCF="T"/>
      </Location>
    </SimpleAllele>
    <RCVList>
      <RCVAccession Accession="RCV000112011" Version="3">
        <RCVClassifications>
          <GermlineClassification>
            <ReviewStatus>criteria provided, single submitter</ReviewStatus>
            <Description DateLastEvaluated="2019-01-01">Likely pathogenic</Description>
          </GermlineClassification>
        </RCVClassifications>
      </RCVAccession>
      <RCVAccession Accession="RCV000031204" Version="5">
        <RCVClassifications>
          <GermlineClassification>
            <ReviewStatus>reviewed by expert panel</ReviewStatus>
            <Description DateLastEvaluated="2016-12-15">Pathogenic</Description>
          </GermlineClassification>
        </RCVClassifications>
      </RCVAccession>
    </RCVList>
    <Classifications>
      {germline}
    </Classifications>
    <ClinicalAssertionList>
      <ClinicalAssertion ID="20157" SubmissionDate="2016-12-15" DateLastUpdated="2017-01-01">
        <ClinVarAccession Accession="SCV000282346" Version="1" SubmitterName="Evidence-based Network for the Interpretation of Germline Mutant Alleles (ENIGMA)"/>
        <RecordStatus>current</RecordStatus>
        <Classification DateLastEvaluated="2016-12-15">
          <ReviewStatus>reviewed by expert panel</ReviewStatus>
          <GermlineClassification>Pathogenic</GermlineClassification>
        </Classification>
        <ObservedInList>
          <ObservedIn>
            <Sample><Origin>germline</Origin></Sample>
            <Method><MethodType>curation</MethodType></Method>
          </ObservedIn>
        </ObservedInList>
      </ClinicalAssertion>
    </ClinicalAssertionList>
  </ClassifiedRecord>
</VariationArchive>
"""

EXPERT_PANEL_GERMLINE = """
      <GermlineClassification DateLastEvaluated="2016-12-15" NumberOfSubmissions="2" NumberOfSubmitters="2" DateCreated="2013-05-01" MostRecentSubmission="2024-02-20">
        <DescriptionHistory Dated="2015-01-01">
          <Description>Likely pathogenic</Description>
        </DescriptionHistory>
        <ReviewStatus>reviewed by expert panel</ReviewStatus>
        <Description>Pathogenic</Description>
        <ConditionList>
          <TraitSet ID="7920" Type="Disease" ContributesToAggregateClassification="true">
            <Trait ID="16761" Type="Disease">
              <Name>
                <ElementValue Type="Alternate">Hereditary breast and ovarian cancer syndrome</ElementValue>
                <XRef ID="D061325" DB="MeSH"/>
              </Name>
              <Name>
                <ElementValue Type="Preferred">Hereditary breast ovarian cancer syndrome</ElementValue>
                <XRef ID="MONDO:0003582" DB="MONDO"/>
              </Name>
              <AttributeSet>
                <Attribute Type="GARD id">15010</Attribute>
                <XRef ID="15010" DB="Office of Rare Diseases"/>
              </AttributeSet>
              <XRef ID="GTR000514601" DB="Genetic Testing Registry (GTR)"/>
              <XRef ID="145" DB="Orphanet"/>
              <XRef ID="C0677776" DB="MedGen"/>
            </Trait>
          </TraitSet>
        </ConditionList>
      </GermlineClassification>
"""

CONFLICTING_GERMLINE = """
      <GermlineClassification DateLastEvaluated="2019-07-02" NumberOfSubmissions="3" NumberOfSubmitters="3" DateCreated="2017-12-26" MostRecentSubmission="2020-06-22">
        <ReviewStatus>criteria provided, conflicting classifications</ReviewStatus>
        <Description>Conflicting classifications of pathogenicity</Description>
        <Explanation DataSource="ClinVar" Type="public">Uncertain significance(2); Likely benign(1)</Explanation>
        <ConditionList>
          <TraitSet ID="7920" Type="Disease" ContributesToAggregateClassification="true">
            <Trait ID="16761" Type="Disease">
              <Name><ElementValue Type="Preferred">Hereditary breast ovarian cancer syndrome</ElementValue></Name>
              <XRef ID="C0677776" DB="MedGen"/>
            </Trait>
          </TraitSet>
          <TraitSet ID="9460" Type="Finding" ContributesToAggregateClassification="true">
            <Trait ID="17556" Type="Finding">
              <Name><ElementValue Type="Preferred">not provided</ElementValue></Name>
              <XRef ID="C3661900" DB="MedGen"/>
            </Trait>
          </TraitSet>
          <TraitSet ID="1234" Type="Disease" ContributesToAggregateClassification="false">
            <Trait ID="999" Type="Disease">
              <Name><ElementValue Type="Preferred">Fanconi anemia</ElementValue></Name>
            </Trait>
          </TraitSet>
        </ConditionList>
      </GermlineClassification>
"""

NO_GERMLINE = """
      <SomaticClinicalImpact>
        <ReviewStatus>no classification provided</ReviewStatus>
      </SomaticClinicalImpact>
"""


def _variation_archive_element(germline):
    return ET.fromstring(VARIATION_ARCHIVE_TEMPLATE.format(germline=germline))


def test_aggregate_classification_single_condition():
    agg = clinvar_common.aggregateClassification(_variation_archive_element(EXPERT_PANEL_GERMLINE))

    assert agg.valid
    assert agg.variationID == '37394'
    assert agg.accession == 'VCV000037394'
    assert agg.version == '12'
    assert agg.dateLastUpdated == '2025-03-01'
    # The current Description, not the one under DescriptionHistory
    assert agg.clinicalSignificance == 'Pathogenic'
    assert agg.reviewStatus == 'reviewed by expert panel'
    assert agg.dateLastEvaluated == '2016-12-15'
    assert agg.numberOfSubmissions == '2'
    assert agg.numberOfSubmitters == '2'
    assert agg.mostRecentSubmission == '2024-02-20'
    assert agg.explanation is None
    assert agg.conditions == ['Hereditary breast ovarian cancer syndrome']
    # Trait and preferred-Name XRefs only: no GTR, no alternate-name MeSH,
    # no AttributeSet Office of Rare Diseases
    assert agg.conditionDbIds == [['Orphanet_145', 'MedGen_C0677776', 'MONDO_MONDO:0003582']]


def test_aggregate_classification_conflicting():
    agg = clinvar_common.aggregateClassification(_variation_archive_element(CONFLICTING_GERMLINE))

    assert agg.clinicalSignificance == 'Conflicting classifications of pathogenicity'
    assert agg.reviewStatus == 'criteria provided, conflicting classifications'
    assert agg.explanation == 'Uncertain significance(2); Likely benign(1)'
    # The TraitSet that doesn't contribute to the aggregate is left out
    assert agg.conditions == ['Hereditary breast ovarian cancer syndrome', 'not provided']
    assert agg.conditionDbIds == [['MedGen_C0677776'], ['MedGen_C3661900']]


def test_aggregate_classification_absent():
    agg = clinvar_common.aggregateClassification(_variation_archive_element(NO_GERMLINE))

    assert not agg.valid
    assert agg.accession == 'VCV000037394'
    assert agg.clinicalSignificance is None
    assert agg.conditions == []


def test_variation_archive_aggregate_is_not_first_rcv():
    va = clinvar_common.variationArchive(_variation_archive_element(EXPERT_PANEL_GERMLINE))

    assert va.valid
    # referenceAssertion takes whichever RCV comes first...
    assert va.referenceAssertion.clinicalSignificance == 'Likely pathogenic'
    # ...while the aggregate is the VCV-level call
    assert va.aggregateClassification.clinicalSignificance == 'Pathogenic'
    assert va.aggregateClassification.reviewStatus == 'reviewed by expert panel'


def _parsed_rows(capsys, germline):
    from . import clinVarParse

    va = clinvar_common.variationArchive(_variation_archive_element(germline))
    clinVarParse.printHeader()
    clinVarParse.processSubmission(va, 'GRCh38')
    lines = capsys.readouterr().out.rstrip('\n').split('\n')
    header = lines[0].split('\t')
    rows = [line.split('\t') for line in lines[1:]]
    for row in rows:
        assert len(row) == len(header)
    return [dict(zip(header, row)) for row in rows]


def test_clinvarparse_writes_aggregate_columns(capsys):
    [row] = _parsed_rows(capsys, CONFLICTING_GERMLINE)

    assert row['SCV'] == 'SCV000282346'
    assert row['ClinicalSignificance'] == 'Pathogenic'
    assert row['VCV_VariationID'] == '37394'
    assert row['VCV_Accession'] == 'VCV000037394'
    assert row['VCV_Version'] == '12'
    assert row['VCV_ClinicalSignificance'] == 'Conflicting classifications of pathogenicity'
    assert row['VCV_ReviewStatus'] == 'criteria provided, conflicting classifications'
    assert row['VCV_DateLastEvaluated'] == '2019-07-02'
    assert row['VCV_NumberOfSubmissions'] == '3'
    assert row['VCV_NumberOfSubmitters'] == '3'
    assert row['VCV_MostRecentSubmission'] == '2020-06-22'
    # ';' would be turned into '.' by convert_tsv_to_vcf.py
    assert row['VCV_Explanation'] == 'Uncertain significance(2), Likely benign(1)'
    assert row['VCV_Conditions'] == 'Hereditary breast ovarian cancer syndrome|not provided'
    assert row['VCV_ConditionDB_IDs'] == 'MedGen_C0677776|MedGen_C3661900'
    assert row['VCV_DateLastUpdated'] == '2025-03-01'


def test_clinvarparse_condition_db_ids_stay_aligned(capsys):
    # The first condition has no IDs; its empty slot must be kept so the
    # second condition's IDs don't shift onto the first
    germline = CONFLICTING_GERMLINE.replace('<XRef ID="C0677776" DB="MedGen"/>', '')
    [row] = _parsed_rows(capsys, germline)

    assert row['VCV_Conditions'] == 'Hereditary breast ovarian cancer syndrome|not provided'
    assert row['VCV_ConditionDB_IDs'] == '|MedGen_C3661900'


def test_clinvarparse_aggregate_columns_absent(capsys):
    [row] = _parsed_rows(capsys, NO_GERMLINE)

    assert row['VCV_VariationID'] == '37394'
    assert row['VCV_Accession'] == 'VCV000037394'
    for column in ('VCV_ClinicalSignificance', 'VCV_ReviewStatus',
                   'VCV_DateLastEvaluated', 'VCV_Explanation',
                   'VCV_Conditions', 'VCV_ConditionDB_IDs'):
        assert row[column] == '-'







