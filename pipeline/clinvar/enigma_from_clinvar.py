import re
from multiprocessing.pool import ThreadPool

import click
import pandas as pd
from lxml import etree

import clinvar
import hgvs_utils

import logging


default_val = 'NaN'  # TODO: just for recreation TODO: what to choose?

MULTI_ENTRY_SEP = ','  # TODO will this break stuff?


def _get_clinvar_sets(fin):
    f = etree.parse(fin)
    root = f.getroot()
    return root.xpath('//ClinVarSet')


def _is_bic_designation(s):
    return any(k in s.lower() for k in {'ins', 'del', 'dup'}) or \
           ('>' in s and not s.startswith('p.')) or \
           s.startswith('p.')


def _fetch_bic(enigma_el):
    elems = enigma_el.xpath(
        "MeasureSet/Measure/Name/ElementValue[@Type='Alternate']")
    elems_text = [e.text.replace(' ', '') for e in elems]

    bic = [e for e in elems_text if _is_bic_designation(e)]

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


def _extract_genomic_coordinates(cvs_el, assembly):
    sequence_location_el = cvs_el.find(
        'ReferenceClinVarAssertion/MeasureSet/Measure/SequenceLocation[@Assembly="{}"]'.format(
            assembly))

    if sequence_location_el is not None:
        return "chr{}:{}:{}>{}".format(
            sequence_location_el.get('Chr'),
            sequence_location_el.get('positionVCF'),
            sequence_location_el.get('referenceAlleleVCF'),
            sequence_location_el.get('alternateAlleleVCF'))

    return default_val


def _extract_condition_info(cvs_el):
    symbol_lst = cvs_el.xpath(
        'ReferenceClinVarAssertion/TraitSet/Trait/Symbol[ElementValue/@Type="Preferred"]')

    if symbol_lst:
        el = symbol_lst[0]

        condition_id_symbol = clinvar.textIfPresent(el, 'ElementValue') # e.g. BROVCA2
        condition_id_type = default_val
        condition_id_id = default_val

        xref_el = el.find('XRef')
        if xref_el is not None:
            condition_id_type = xref_el.get('DB') # e.g. OMIM
            condition_id_id = xref_el.get('ID') # e.g. 612555

    omim_name_el = cvs_el.xpath('ReferenceClinVarAssertion/TraitSet/Trait/Name[XRef/@DB="OMIM"]')

    if omim_name_el:
        id_value = clinvar.textIfPresent(omim_name_el[0], "ElementValue")
        condition_id_value = "{}; {} ({})".format(id_value, condition_id_symbol, condition_id_id)
    else:
        condition_id_type = default_val

    return condition_id_type, condition_id_value


def _xpath(el, xpath):
    e = el.xpath(xpath)
    if e:
        return e[0]
    return default_val

# clinvar set element
def parse_record(cvs_el, hgvs_util, assembly="GRCh38"):
    rec = {}

    # TODO: careful different logic. processing happens on variant set instead per assertion as in main clinvar parsing. make it more explicits?
    enigma_assertion = cvs_el.xpath('ClinVarAssertion[contains(ClinVarSubmissionID/@submitter, "ENIGMA")]')[0]

    rec["Gene_symbol"] = clinvar.textIfPresent(cvs_el,
                                               'ReferenceClinVarAssertion/MeasureSet/Measure/MeasureRelationship/Symbol/ElementValue[@Type="Preferred"]')

    rec["Genomic_Coordinate"] = _extract_genomic_coordinates(cvs_el, assembly)

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

    rec["BIC_Nomenclature"] = _fetch_bic(enigma_assertion)

    rec["Abbrev_AA_change"], rec["HGVS_protein"] = _compute_protein_changes(
        hgvs_cdna_complete, hgvs_util)

    rec["URL"] = clinvar.textIfPresent(enigma_assertion,
                                       'ClinicalSignificance/Citation/URL')

    rec["Condition_ID_type"], rec["Condition_ID_value"] = _extract_condition_info(cvs_el)

    trait_set_el = cvs_el.find('ReferenceClinVarAssertion/TraitSet')
    if trait_set_el is not None:
        rec["Condition_category"] = trait_set_el.get('Type')
    else:
        rec["Condition_category"] = default_val

    rec["Clinical_significance"] = clinvar.textIfPresent(enigma_assertion,
                                                         "ClinicalSignificance/Description")

    rec["Date_last_evaluated"] = _xpath(enigma_assertion, "ClinicalSignificance/@DateLastEvaluated")

    # TODO: prefix "PMID" ? handle multiple case
    # TODO: handle special case if empty?
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

    rec["Allele_origin"] = clinvar.textIfPresent(enigma_assertion,
                                                 "ObservedIn/Sample/Origin")

    clinvar_accession_el = enigma_assertion.find("ClinVarAccession")
    # TODO include accession version rec["ClinVarAccession"] = "{}.{}".format(clinvar_accession_el.get("Acc"), clinvar_accession_el.get("Version"))

    rec["ClinVarAccession"] = clinvar_accession_el.get("Acc")

    return rec


@click.command()
@click.argument('filtered_clinvar_xml', type=click.Path(exists=True))
@click.argument('output', type=click.Path())
def main(filtered_clinvar_xml, output):
    enigma_sets = _get_clinvar_sets(filtered_clinvar_xml)

    # TODO: tweak, necessary?
    logging.basicConfig(level=logging.INFO)

    # TODO: reshuffle to use multiprocessing?
    # TODO: make sure stack trace gets up in case of errors
    # from multiprocessing import Pool
    # sets_str = [etree.tostring(el, encoding='unicode', pretty_print=False) for el in sets]
    # pool = Pool(6)
    # variant_records = pool.map(parse_record_str, sets_str)

    pool = ThreadPool(6)
    # variant_records = [parse_record_str(s) for s in sets_str[0:100]]

    # TODO: fix naming?
    hgvs_util = hgvs_utils.HGVSWrapper()

    variant_records = pool.map(lambda s: parse_record(s, hgvs_util),
                               enigma_sets)

    df = pd.DataFrame.from_dict(variant_records)

    print(df.info())

    df.to_csv(output, sep='\t', index=False)

    # TODO: set proper return code!


if __name__ == "__main__":
    main()
