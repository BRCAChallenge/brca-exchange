import copy
import itertools
import logging
import re
import click
import pandas as pd
import xml.etree.ElementTree as ET


import common
from clinvar import clinvar_common as clinvar
from clinvar import hgvs_utils

default_val = None

MULTI_ENTRY_SEP = ','


def _get_variation_archives(fin):
    f = etree.parse(fin)
    root = f.getroot()
    return root.xpath('//VariationArchive')


three_letters_aa = re.compile('p.\\(?[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}') # e.g. p.(Tyr831SerfsTer9)
def _is_bic_designation(s):
    return any(k in s.lower() for k in {'ins', 'del', 'dup'}) or \
        (not s.startswith('p.') and '>' in s) or \
        (s.startswith('p.') and ':' not in s and three_letters_aa.match(s) is None) # shouldn't match for a BIC designator


#
# 11/27/2024: ClinVar no longer supplies the BIC nomenclature terms in their records
def _fetch_bic(cvs_el):
    return default_val


def _compute_protein_changes(hgvs_cdna, hgvs_util):
    if hgvs_cdna is not None:
        v_protein = hgvs_util.compute_protein_change(hgvs_cdna)
        if v_protein is not None:
            hgvs_protein = v_protein.format(
                {'p_3_letter': True, 'p_term_asterisk': False})
            abbrev_aa_change = v_protein.format(
                {'p_3_letter': False, 'p_term_asterisk': True})
            if ':' in hgvs_protein:
                hgvs_protein = hgvs_protein.split(":")[1]
            else:
                hgvs_protein = default_val
            if ':' in abbrev_aa_change:
                abbrev_aa_change = abbrev_aa_change.split(":")[1]. \
                    lstrip("p.").replace('(', '').replace(')', '')
            else:
                abbrev_aa_change = default_val
            return (abbrev_aa_change, hgvs_protein)
    return (default_val, default_val)





def _xpath(el, xpath):
    e = el.xpath(xpath)
    if e:
        return e[0]
    return default_val


def _xpath_text(el, xpath):
    e = _xpath(el, xpath)
    if e is not None:
        return e.text
    return default_val


def _parse_engima_assertion(enigma_assertion, hgvs_util):
    rec = {}

    rec["URL"] = default_val
    rec["Clinical_significance"] = enigma_assertion.clinicalSignificance
    rec["Date_last_evaluated"] = enigma_assertion.dateSignificanceLastEvaluated
    rec["Assertion_method"] = enigma_assertion.assertionMethod
    rec["Assertion_method_citation"] = enigma_assertion.assertionMethodCitation
    rec["Comment_on_clinical_significance"] = ','.join(enigma_assertion.description)
    rec["Collection_method"] = ','.join(enigma_assertion.method)
    rec["Allele_origin"] = ','.join(enigma_assertion.origin)
    rec["ClinVarAccession"] = "%s.%s" (enigma_assertion.accession, enigma_assertion.accession_version)

    return rec

# Parse a ClinVar variation archive element
def parse_record(va_el, hgvs_util, symbols, assembly="GRCh38"):
    '''
    Extracts information out of a VariationArchive element

    :param cvs_el: VariationArchive element
    :param hgvs_util: instance of HGVS util to calculate protein changes
    :param assembly: str assembly to use
    :return: list of dictionary containing the extracted data. keys correspond to column names.
    Each element of the list corresponds to a ENIGMA submission
    '''
    rec = {}
    va = clinvar.variationArchive(va_el)
    variant = va.variant
    
    rec["Gene_symbol"] = variant.geneSymbol
    if assembly not in variant.coordinates.keys():
        logging.warning("Skipping variant %s as no genomic coordinates could be extracted", va.name)
        return []
    coords = variant.coordinates[assembly]
    rec["Genomic_Coordinate"] = "chr%d:%d:%s>%s" % (coords.chr, coords.pos, coords.ref, coords.alt)
    rec["BIC_Nomenclature"] = _fetch_bic(va_el)
    rec["Reference_sequence"] = default_val
    rec["HGVS_cDNA"] = default_val
    rec["Abbrev_AA_change"] = default_val
    rec["HGVS_protein"] = default_val
    if re.match("NM_", variant.hgvs_cdna):
        rec["Reference_sequence"] = variant.hgvs_cdna.split(":")[0]
        rec["HGVS_cDNA"] = variant.hgvs_cdna.split(":")[1]
        (rec["Abbrev_AA_change"], rec["HGVS_protein"]) = _compute_protein_changes(variant.hgvs_cdna, hgvs_util)
    rec["Condition_ID_type"] = va.classification.condition_type
    rec["Condition_ID_value"] = va.classification.condition_value
    #
    # 11/25/2024: the trait set / condition category has been deprecated (is no longer displayed)
    rec["Condition_category"] = default_val
    for scv_accession in va.otherAssertions.keys():
        oa = va.otherAssertions[scv_accession]
        if oa.reviewStatus == "reviewed by expert panel":
        rec_tmp = copy.deepcopy(rec)
        rec_tmp.update(_parse_engima_assertion(oa, hgvs_util))
        return(rec_tmp)
    return None


def _create_df(variant_records):
    df = pd.DataFrame.from_dict(variant_records)

    df['BX_ID'] = pd.Series(range(1, len(df) + 1))

    target_header = ['Gene_symbol',
                     'Genomic_Coordinate',
                     'Reference_sequence',
                     'HGVS_cDNA',
                     'BIC_Nomenclature',
                     'Abbrev_AA_change',
                     'URL',
                     'Condition_ID_type',
                     'Condition_ID_value',
                     'Condition_category',
                     'Clinical_significance',
                     'Date_last_evaluated',
                     'Assertion_method',
                     'Assertion_method_citation',
                     'Clinical_significance_citations',
                     'Comment_on_clinical_significance',
                     'Collection_method',
                     'Allele_origin',
                     'ClinVarAccession',
                     'HGVS_protein',
                     'BX_ID']

    if df.empty:
        # can happen if none of the genes of interest are curated by ENIGMA
        return pd.DataFrame({}, columns=target_header)

    return df.loc[:, target_header]


@click.command()
@click.argument('filtered_clinvar_xml', type=click.Path(exists=True))
@click.argument('output', type=click.Path(writable=True))
@click.option('--logs', type=click.Path(writable=True))
@click.option('--gene', type=str, required=True, multiple=True)
def main(filtered_clinvar_xml, gene, output, logs):
    common.utils.setup_logfile(logs)
    hgvs_util = hgvs_utils.HGVSWrapper()
    gene_symbols = list(set(gene))
    variant_records = list()
    with open(filtered_clinvar_xml) as inputFile:
        for event, elem in ET.iterparse(inputFile, events=('start', 'end')):
            if event == 'end' and elem.tag == 'VariationArchive':
                next_record = parse_record(elem, hgvs_util, gene_symbol)
                variant_records.append(next_record)
    df = _create_df(variant_records)
    df.to_csv(output, sep='\t', index=False, na_rep='-')


if __name__ == "__main__":
    main()
