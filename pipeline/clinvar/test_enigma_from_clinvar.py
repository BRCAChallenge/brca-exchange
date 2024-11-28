import os

import pytest
from lxml import etree

from . import enigma_from_clinvar
from . import hgvs_utils


@pytest.mark.parametrize("s,expected", [
    ['2415delAG', True],
    ['2829G>T', True],
    ['p.Gln494*', True],
    ['p.(Tyr831SerfsTer9)', False],
    ['p.K862E:AAG>GAG', False],
    ['Q1313X', False]
])

def test_parse_record(hgvs_wrapper):
    hgvs_util = hgvs_utils.HGVSWrapper(hgvs_wrapper.hgvs_dp)

    xml_path = os.path.join(os.path.dirname(__file__), 'test_files',
                            'enigma_clinvar_set.xml')
    cvs_el = etree.parse(xml_path)

    res = enigma_from_clinvar.parse_record(cvs_el, hgvs_util, ['BRCA1', 'BRCA2'])

    expected = [{
        'Allele_origin': 'Germline',
        'Assertion_method': 'ENIGMA BRCA1/2 Classification Criteria (2017-06-29)',
        'Assertion_method_citation': 'https://submit.ncbi.nlm.nih.gov/ft/byid/vuhooppz/enigma_rules_2017-06-29-v2.5.1.pdf',
        'BIC_Nomenclature': None,
        'ClinVarAccession': 'SCV000783591.1',
        'Clinical_significance': 'Pathogenic',
        'Clinical_significance_citations': '',
        'Collection_method': 'Curation',
        'Comment_on_clinical_significance': 'Variant allele predicted to encode a truncated non-functional protein.',
        'Date_last_evaluated': '2018-07-13',
        'Gene_symbol': 'BRCA1',
        'Genomic_Coordinate': 'chr17:43094316:TTGATTCAGACTCCCCATCATGTGAGTCATCAGAACCTAACA>T',
        'HGVS_cDNA': 'c.1175_1215del',
        'HGVS_protein': 'p.Leu392fs',
        'Reference_sequence': 'NM_007294.4',
        'URL': None,
    }]

    assert res == expected
