import os

import pytest
from lxml import etree

import enigma_from_clinvar
import hgvs_utils


@pytest.mark.parametrize("s,expected", [
    ['2415delAG', True],
    ['2829G>T', True],
    ['p.Gln494*', True],
    ['p.(Tyr831SerfsTer9)', False],
    ['p.K862E:AAG>GAG', False],
    ['Q1313X', False]
])
def test_is_bic_designation(s, expected):
    assert enigma_from_clinvar._is_bic_designation(s) == expected


def test_parse_record():
    hgvs_util = hgvs_utils.HGVSWrapper()

    xml_path = os.path.join(os.path.dirname(__file__), 'test_files',
                            'enigma_clinvar_set.xml')
    cvs_el = etree.parse(xml_path)

    res = enigma_from_clinvar.parse_record(cvs_el, hgvs_util)

    expected = [{
        'Abbrev_AA_change': u'I562Mfs*12',
        'Allele_origin': 'Germline',
        'Assertion_method': 'ENIGMA BRCA1/2 Classification Criteria (2017-06-29)',
        'Assertion_method_citation': 'https://submit.ncbi.nlm.nih.gov/ft/byid/ncs3jnil/enigma_rules_2017-06-29.pdf',
        'BIC_Nomenclature': None,
        'ClinVarAccession': 'SCV000783130.1',
        'Clinical_significance': 'Pathogenic',
        'Clinical_significance_citations': '',
        'Collection_method': 'Curation',
        'Comment_on_clinical_significance': 'Variant allele predicted to encode a truncated non-functional protein.',
        'Condition_ID_type': 'OMIM',
        'Condition_ID_value': 'BREAST-OVARIAN CANCER, FAMILIAL, SUSCEPTIBILITY TO, 1; BROVCA1 (604370)',
        'Condition_category': 'Disease',
        'Date_last_evaluated': '2017-12-15',
        'Gene_symbol': 'BRCA1',
        'Genomic_Coordinate': 'chr17:43093845:A>ACTTTC',
        'HGVS_cDNA': 'c.1685_1686insGAAAG',
        'HGVS_protein': u'p.(Ile562MetfsTer12)',
        'Reference_sequence': 'NM_007294.3',
        'URL': None,
    }]

    assert res == expected
