#!/usr/bin/env python
"""
enigma_postprocess.py: postprocess the ENIGMA data to fix data items such as
broken URLs or variants with specific issues
"""
import argparse
import csv
import logging

def classificationToBeRemoved(row):
    """
    The variant NM_000059.3:c.476-2A>G has classifications by ENIGMA as both 
    pathogenic and uncertain.  The ENIGMA consortium wants to retain only the
    pathogenic classification.  Return a binary flag indicating if this is
    that classification, and should be skipped..
    """
    if row["Reference_sequence"] == "NM_000059.3" \
            and row["HGVS_cDNA"] == "c.476-2A>G" \
            and row["Clinical_significance"] != "Pathogenic":
        return True
    else:
        return False


def fillInFounderMutations(row):
    """
    Make sure that the founder mutations are represented with each of their
    BIC nomenlcature terms.
    The founder mutations are
    BRCA1 185_186delAG (a.k.a. 185delAG, 186delAG, 187delAG)
    BRCA1 5382_5383insC (a.k.a. 5382insC, 5383insC, 5384insC, 5385insC)
    BRCA2 6174delT
    Here, the first and second columns represent the gene name and BIC 
    nomenclature terms given by ENIGMA, and in parentheses are the other
    BIC terms to be added.  The single BIC term can be turned into a list of
    BIC terms by delimiting the list with commas.
    Addresses GitHub Issue #213 
    (https://github.com/BD2KGenomics/brca-exchange/issues/213)
    """
    if row["Gene_symbol"] == "BRCA1":
        if row["BIC_Nomenclature"] == "185_186delAG":
            row["BIC_Nomenclature"] += ",185delAG|186delAG,187delAG"
        elif row["BIC_Nomenclature"] == "5382_5383insC":
            row["BIC_Nomenclature"] += "|5382insC,5383insC,5384insC,5385insC"
    return row


def fixAssertionCitation(row):
    """
    The Assertion Citation refers to the set of rules that were used by
    ENIGMA in their classification.  These rules documents are not necessarily
    accessible from outside the ENIGMA site, i.e. you can't link to them 
    directly from an outside URL, but if you link directly to a general URL,
    you can navigate from there to the rules document.
    Addresses GitHub Issue #237
    (https://github.com/BD2KGenomics/brca-exchange/issues/237)
    """
    row["Assertion_method_citation"] = "https://enigmaconsortium.org/library/general-documents/"
    return row


def fixBrokenAminoAcidChange(row):
    """
    For one particular variant, NM_000059.3(BRCA1):p.(Arg1190Trp), the one-letter
    amino acid change is listed by ENIGMA as P1190S but should be R1190W.  A user noticed
    the issue.  Mandy Spurdle confirms that this should be changed.
    Addresses GitHub Issue # 184
    (https://github.com/BD2KGenomics/brca-exchange/issues/184)
    """
    if row["HGVS_protein"] == "p.(Arg1190Trp)" and row["Abbrev_AA_change"] == "P1190S":
        row["Abbrev_AA_change"] = "R1190W"
    return row

def updateNotYetClassifiedVariants(row):
    """
    For variants with the expert-reviewed clinical status of "Not Yet Classified",
    change the status to "Not Yet Reviewed"
    Addresses GitHub Issue # 362
    (https://github.com/BD2KGenomics/brca-exchange/issues/184)
    """
    if row["Clinical_significance"] == "Not Yet Classified":
        row["Clinical_significance"] = "Not Yet Reviewed"
    return row

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", 
                        help="input file with merged ENIGMA data")
    parser.add_argument("-o", "--output",
                        help="Output file with corrected ENIGMA data")
    parser.add_argument('-a', "--artifacts_dir", 
                        help='Artifacts directory with pipeline artifact files.')
    parser.add_argument("-v", "--verbose", action="count", default=False, 
                        help="determines logging")
    args = parser.parse_args()

    if args.verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.CRITICAL
    log_file_path = args.artifacts_dir + "enigma-postprocess.log"
    logging.basicConfig(filename=log_file_path, filemode="w", level=logging_level)

    csvIn = csv.DictReader(open(args.input, "r"), delimiter='\t')
    csvOut = csv.DictWriter(open(args.output, "w"), delimiter='\t',
                            fieldnames=csvIn.fieldnames)
    csvOut.writerow(dict((fn, fn) for fn in csvIn.fieldnames))
    for row in csvIn:
        row = fillInFounderMutations(row)
        row = fixAssertionCitation(row)
        row = fixBrokenAminoAcidChange(row)
        if not classificationToBeRemoved(row):
            csvOut.writerow(row)


if __name__ == "__main__":
    main()
