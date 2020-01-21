import copy
import itertools
import re

import click
import pandas as pd
from lxml import etree

import clinvar
import common
import hgvs_utils

default_val = None

MULTI_ENTRY_SEP = ','


def _get_clinvar_sets(fin):
    f = etree.parse(fin)
    root = f.getroot()
    return root.xpath('//ClinVarSet')


three_letters_aa = re.compile('p.\(?[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}') # e.g. p.(Tyr831SerfsTer9)
def _is_bic_designation(s):
    return any(k in s.lower() for k in {'ins', 'del', 'dup'}) or \
        (not s.startswith('p.') and '>' in s) or \
        (s.startswith('p.') and ':' not in s and three_letters_aa.match(s) is None) # shouldn't match for a BIC designator


def _fetch_bic(cvs_el):
    elems = cvs_el.xpath(
        "*/MeasureSet/Measure/Name/ElementValue[@Type='Alternate']")
    elems_text = [e.text for e in elems]

    bic = sorted({e for e in elems_text if _is_bic_designation(e)})

    if bic:
        return '|'.join(bic)
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

            return abbrev_aa_change, hgvs_protein

    return default_val, default_val


def _extract_assertion_method(enigma_assertion):
    method_el = enigma_assertion.find(
        "AttributeSet/Attribute[@Type='AssertionMethod']")

    if method_el is not None:
        attribute_set_el = method_el.getparent()
        return method_el.text, clinvar.textIfPresent(attribute_set_el,
                                                     "Citation/URL")
    else:
        return default_val, default_val


def _extract_condition_info(cvs_el):
    symbol_lst = cvs_el.xpath(
        'ReferenceClinVarAssertion/TraitSet/Trait/Symbol[ElementValue/@Type="Preferred"]')

    condition_id_type = default_val

    if symbol_lst:
        el = symbol_lst[0]

        condition_id_symbol = clinvar.textIfPresent(el,
                                                    'ElementValue')  # e.g. BROVCA2
        condition_id_id = default_val

        xref_el = el.find('XRef')
        if xref_el is not None:
            condition_id_type = xref_el.get('DB')  # e.g. OMIM
            condition_id_id = xref_el.get('ID')  # e.g. 612555

    omim_name_el = cvs_el.xpath(
        'ReferenceClinVarAssertion/TraitSet/Trait/Name[XRef/@DB="OMIM"]')

    if omim_name_el:
        id_value = clinvar.textIfPresent(omim_name_el[0], "ElementValue")
        condition_id_value = "{}; {} ({})".format(id_value, condition_id_symbol,
                                                  condition_id_id)
    else:
        condition_id_value = default_val

    return condition_id_type, condition_id_value


def _extract_clinvar_accession(engima_assertion):
    clinvar_accession_el = engima_assertion.find("ClinVarAccession")

    if clinvar_accession_el is not None:
        acc = clinvar_accession_el.get("Acc")
        v = clinvar_accession_el.get("Version")

        if acc is not None:
            if v is not None:
                return "{}.{}".format(acc, v)
            else:
                return acc

    return default_val


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

    # getting hgvs cdna
    hgvs_list = [e.text for e in enigma_assertion.findall(
        'MeasureSet/Measure/AttributeSet/Attribute[@Type="HGVS"]')]

    hgvs_list_filtered = [s for s in hgvs_list if re.match('NM_.*:.*', s)]

    if hgvs_list_filtered:
        hgvs_cdna_complete = hgvs_list_filtered[0]
        rec["Reference_sequence"] = hgvs_cdna_complete.split(":")[0]
        rec["HGVS_cDNA"] = hgvs_cdna_complete.split(":")[1]
    else:
        hgvs_cdna_complete = None
        rec["Reference_sequence"] = default_val
        rec["HGVS_cDNA"] = default_val

    rec["Abbrev_AA_change"], rec["HGVS_protein"] = _compute_protein_changes(
        hgvs_cdna_complete, hgvs_util)

    rec["URL"] = clinvar.textIfPresent(enigma_assertion,
                                       'ClinicalSignificance/Citation/URL')

    rec["Clinical_significance"] = clinvar.textIfPresent(enigma_assertion,
                                                         "ClinicalSignificance/Description")

    rec["Date_last_evaluated"] = _xpath(enigma_assertion,
                                        "ClinicalSignificance/@DateLastEvaluated")

    rec[
        "Clinical_significance_citations"] = MULTI_ENTRY_SEP.join(
        ["PMID:{}".format(e.text) for e in enigma_assertion.findall(
            "ClinicalSignificance/Citation/ID[@Source='PubMed']")])

    rec["Comment_on_clinical_significance"] = clinvar.textIfPresent(
        enigma_assertion, "ClinicalSignificance/Comment").replace('\t', '    ')

    rec["Assertion_method"], rec[
        "Assertion_method_citation"] = _extract_assertion_method(
        enigma_assertion)

    rec["Collection_method"] = clinvar.textIfPresent(enigma_assertion,
                                                     "ObservedIn/Method/MethodType")
    if rec["Collection_method"] is not None:
        rec["Collection_method"] = rec["Collection_method"].capitalize()

    rec["Allele_origin"] = clinvar.textIfPresent(enigma_assertion,
                                                 "ObservedIn/Sample/Origin")
    if rec["Allele_origin"] is not None:
        rec["Allele_origin"] = rec["Allele_origin"].capitalize()

    rec["ClinVarAccession"] = _extract_clinvar_accession(enigma_assertion)

    return rec

# clinvar set element
def parse_record(cvs_el, hgvs_util, assembly="GRCh38"):
    '''
    Extracts information out of a ClinVarSet XML element

    :param cvs_el: ClinVarSet XML element
    :param hgvs_util: instance of HGVS util to calculate protein changes
    :param assembly: str assembly to use
    :return: list of dictionary containing the extracted data. keys correspond to column names.
    Each element of the list corresponds to a ENIGMA submission
    '''
    rec = {}

    rec["Gene_symbol"] = _xpath_text(cvs_el,
                                     'ReferenceClinVarAssertion/MeasureSet/Measure/MeasureRelationship/Symbol/ElementValue[starts-with(., "BRCA") and @Type="Preferred"]')

    measure_el = cvs_el.find('ReferenceClinVarAssertion/MeasureSet/Measure')
    rec["Genomic_Coordinate"] = str(clinvar.extract_genomic_coordinates_from_measure(measure_el)[assembly]).replace('g.', '')

    rec["BIC_Nomenclature"] = _fetch_bic(cvs_el)

    rec["Condition_ID_type"], rec[
        "Condition_ID_value"] = _extract_condition_info(cvs_el)

    trait_set_el = cvs_el.find('ReferenceClinVarAssertion/TraitSet')
    if trait_set_el is not None:
        rec["Condition_category"] = trait_set_el.get('Type')
    else:
        rec["Condition_category"] = default_val

    lst = []
    for enigma_assertion in cvs_el.xpath(
        'ClinVarAssertion[contains(ClinVarSubmissionID/@submitter, "ENIGMA")]'):

        rec_tmp = copy.deepcopy(rec)
        rec_tmp.update(_parse_engima_assertion(enigma_assertion, hgvs_util))
        lst.append(rec_tmp)

    return lst


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

    return df.loc[:, target_header]


@click.command()
@click.argument('filtered_clinvar_xml', type=click.Path(exists=True))
@click.argument('output', type=click.Path(writable=True))
@click.option('--logs', type=click.Path(writable=True))
def main(filtered_clinvar_xml, output, logs):
    common.utils.setup_logfile(logs)

    enigma_sets = _get_clinvar_sets(filtered_clinvar_xml)

    hgvs_util = hgvs_utils.HGVSWrapper()

    variant_records_lsts = [ parse_record(s, hgvs_util) for s in enigma_sets ]

    # flattening list of lists
    variant_records = list(itertools.chain.from_iterable(variant_records_lsts))
    df = _create_df(variant_records)
    df.to_csv(output, sep='\t', index=False, na_rep='-')


if __name__ == "__main__":
    main()
