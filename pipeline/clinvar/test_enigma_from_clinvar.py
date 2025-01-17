import os

import pytest
import xml.etree.ElementTree as ET
import sys

current_dir = os.getcwd()
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
import common
import common.hgvs_utils



from . import enigma_from_clinvar




def test_parse_record():
    hgvs_util = common.hgvs_utils.HgvsWrapper()

    xml_path = os.path.join(os.path.dirname(__file__), 'test_files',
                            'enigma_clinvar_set.xml')
    tree = ET.parse(open(xml_path))
    root = tree.getroot()
    for element in root:
        res = enigma_from_clinvar.parse_record(element, hgvs_util, ['BRCA1', 'BRCA2'], {'BRCA1':'NM_007294.4', 'BRCA2': 'NM_000059.4'})
        expected = {
            'Abbrev_AA_change': 'L392Qfs*6',
            'Allele_origin': 'germline',
            'Assertion_method': 'ENIGMA BRCA1/2 Classification Criteria (2017-06-29)',
            'Assertion_method_citation': 'https://submit.ncbi.nlm.nih.gov/ft/byid/vuhooppz/enigma_rules_2017-06-29-v2.5.1.pdf',
            'BIC_Nomenclature': None,
            'ClinVarAccession': 'SCV000783591.1',
            'Clinical_significance': 'Pathogenic',
            'Clinical_significance_citations': '',
            'Collection_method': 'curation',
            'Comment_on_clinical_significance': 'Variant allele predicted to encode a truncated non-functional protein.',
            'Condition_category': None,
            'Condition_ID_type': 'Disease',
            'Condition_ID_value': 'Breast-ovarian cancer, familial, susceptibility to, 1',
            'Date_last_evaluated': '2017-12-15',
            'Gene_symbol': 'BRCA1',
            'Genomic_Coordinate': 'chr17:43094315:TTGATTCAGACTCCCCATCATGTGAGTCATCAGAACCTAACA>T',
            'HGVS_cDNA': 'c.1175_1215del',
            'HGVS_protein': 'p.(Leu392GlnfsTer6)',
            'Reference_sequence': 'NM_007294.4',
            'URL': None,
        }
        assert res == expected
