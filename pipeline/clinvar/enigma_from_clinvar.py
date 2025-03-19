import copy
import itertools
import logging
import re
import click
import pandas as pd
import xml.etree.ElementTree as ET


import common
import common.hgvs_utils
import common.ucsc
from clinvar import clinvar_common as clinvar


default_val = None

MULTI_ENTRY_SEP = ','


three_letters_aa = re.compile('p.\\(?[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}') # e.g. p.(Tyr831SerfsTer9)

def _is_bic_designation(s):
    return any(k in s.lower() for k in {'ins', 'del', 'dup'}) or \
        (not s.startswith('p.') and '>' in s) or \
        (s.startswith('p.') and ':' not in s and three_letters_aa.match(s) is None) # shouldn't match for a BIC designator


def _compute_protein_changes(hgvs_cdna, hgvs_util):
    if hgvs_cdna is not None:
        v_protein = hgvs_util.cdna_to_protein(hgvs_cdna, return_str=False) 
        if v_protein is not None:
            hgvs_protein = str(v_protein.format(
                {'p_3_letter': True, 'p_term_asterisk': False}))
            abbrev_aa_change = str(v_protein.format(
                {'p_3_letter': False, 'p_term_asterisk': True}))
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


def _fetch_bic(va_el):
    bic_list = list()
    for name_item in va_el.findall("./ClassifiedRecord/SimpleAllele/OtherNameList/Name"):
        name = name_item.text
        if _is_bic_designation(name):
            bic_list.append(name)
    if len(bic_list) > 0:
        return '|'.join(bic_list)
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
    rec["URL"] = default_val
    rec["Clinical_significance"] = enigma_assertion.clinicalSignificance
    rec["Clinical_significance_citations"] = ','.join(enigma_assertion.clinicalSignificanceCitations)
    rec["Date_last_evaluated"] = enigma_assertion.dateSignificanceLastEvaluated
    rec["Assertion_method"] = enigma_assertion.assertionMethod
    rec["Assertion_method_citation"] = enigma_assertion.assertionMethodCitation
    rec["Comment_on_clinical_significance"] = enigma_assertion.summaryEvidence
    rec["Collection_method"] = ','.join(enigma_assertion.method).capitalize()
    rec["Allele_origin"] = ','.join(enigma_assertion.origin).capitalize()
    rec["ClinVarAccession"] = "%s.%s" % (enigma_assertion.accession,
                                         enigma_assertion.accession_version)
    return rec

# Parse a ClinVar variation archive element
def parse_record(va_el, hgvs_util, symbols, mane_transcript,
                 assembly="GRCh38", debug=False):
    '''
    Extracts information out of a VariationArchive element
    :param va_el: VariationArchive element
    :param hgvs_util: instance of HGVS util to calculate protein changes
    :param symbols: gene symbols to be parsed for
    :param mane_transcripts: dict mapping symbols to MANE Select transcripts
    :param assembly: str assembly to use
    :return: list of dictionary containing the extracted data. keys correspond to column names.
    Each element of the list corresponds to a ENIGMA submission
    '''
    rec = {}
    va = clinvar.variationArchive(va_el, gene_symbols=symbols, debug=debug,
                                  mane_transcripts=mane_transcript)
    if not hasattr(va, 'variant'):
        logging.warning("Skipping VariationArchive %s as no variant record was found", va.id)
        return(None)
    variant = va.variant
    if variant.valid:
        rec["Gene_symbol"] = variant.geneSymbol
        if assembly not in variant.coordinates.keys():
            logging.warning("Skipping variant %s as no genomic coordinates could be extracted", va.name)
            return None
        coords = variant.coordinates[assembly]
        rec["Genomic_Coordinate"] = "chr%d:%d:%s>%s" % (coords.chr, coords.pos, coords.ref, coords.alt)
        rec["BIC_Nomenclature"] = _fetch_bic(va_el)
        rec["Reference_sequence"] = default_val
        rec["HGVS_cDNA"] = default_val
        rec["Abbrev_AA_change"] = default_val
        rec["HGVS_protein"] = default_val
        #
        # Get the MANE Select cDNA transcript for the reference sequence
        # and HGVS cDNA and protein.
        # 1. If it's already the variant name, do nothing (common case)
        # 2. Else, if it's any of the variant synonyms, use that
        # 3. Else, translate the named transcript to genomic HGVS, and
        #    from that to the MANE Select transcript.
        transcript = None
        hgvs_obj = None
        target_reference_sequence = mane_transcript[rec["Gene_symbol"]]
        if re.search("^" + target_reference_sequence, va.name):
            transcript = va.name
        else:
            if debug:
                print("Looking for transcript in synonyms")
            for synonym in variant.synonyms:
                if re.search("^" + target_reference_sequence, synonym):
                    transcript = synonym
                    break
            if transcript is None:
                if debug:
                    print("Trying to translate HGVS string", va.name)
                hgvs_obj = hgvs_util.parse_hgvs_string(va.name)
                if not hgvs_obj:
                    logging.warning("Skipping variant %s because its HGVS could nt be parsed", va.name)
                    return(None)
                transcript_obj = hgvs_util.to_cdna(hgvs_obj,
                                                   target_transcript=target_reference_sequence)
                if not transcript_obj:
                    logging.warning("Skipping variant %s because its HGVS could not be mapped", va.name)
                    return(None)
                transcript = str(transcript_obj)
        if transcript:
            if not hgvs_obj:
                hgvs_obj = hgvs_util.parse_hgvs_string(transcript)
                if not hgvs_obj:
                    logging.warning("Skipping variant %s due to HGVS parsing error", va.name)
                    return(None)
            (rec["Abbrev_AA_change"], rec["HGVS_protein"]) = _compute_protein_changes(hgvs_obj, hgvs_util)
            rec["Reference_sequence"] = transcript.split(":")[0]
            rec["HGVS_cDNA"] = transcript.split(":")[1]
        rec["Condition_ID_type"] = va.classification.condition_type
        if va.classification.condition_value == None:
            rec["Condition_ID_value"] = "not provided"
        else:
            rec["Condition_ID_value"] = va.classification.condition_value
        #
        # 2/28/25: it seems counter-intuitive that the condition category
        # comes from the ClinVar condition type field, while the condition
        # type is blank, but that's the way it is.  FWIW, neither of these
        # fields are actually displayed.
        rec["Condition_category"] = va.classification.condition_type
        rec["Condition_type"] = default_val
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
    hgvs_util = common.hgvs_utils.HgvsWrapper()
    gene_symbols = list(set(gene))
    mane_transcripts = common.ucsc.symbols_to_mane_transcripts(gene_symbols)
    variant_records = list()
    with open(filtered_clinvar_xml) as inputFile:
        for event, elem in ET.iterparse(inputFile, events=('start', 'end')):
            if event == 'end' and elem.tag == 'VariationArchive':
                next_record = parse_record(elem, hgvs_util, gene_symbols,
                                           mane_transcripts, debug=False)
                if next_record is not None:
                    variant_records.append(next_record)
                elem.clear()
    df = _create_df(variant_records)
    df.to_csv(output, sep='\t', index=False, na_rep='-')


if __name__ == "__main__":
    main()
